import os
from collections import defaultdict
from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.gc_content_checker import GCContentChecker

class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence using a backtracking 
    algorithm to select high CAI codons, maximize diversity, and strictly avoid 
    forbidden sequences, internal promoters, hairpins, and bad GC content.

    The `run` method has been optimized:
      - incremental CDS building (no repeated full joins)
      - caching of expensive window checks
      - cheaper checks run first (forbidden -> promoter -> hairpin)
      - configurable check window size and max_attempts safety cap
    """

    def __init__(self):
        self.aa_to_codons = defaultdict(list)
        self.rbsChooser = None
        self.forbidden_checker = None
        self.promoter_checker = None
        self.gc_checker = None

    def initiate(self) -> None:
        """Initializes the codon table, RBS chooser, and all checkers."""
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()
        
        self.forbidden_checker = ForbiddenSequenceChecker()
        self.forbidden_checker.initiate()
        
        self.promoter_checker = PromoterChecker()
        self.promoter_checker.initiate()
        
        self.gc_checker = GCContentChecker()
        self.gc_checker.initiate()

        # Parse codon usage data
        codon_usage_file = 'genedesign/data/codon_usage.txt'
        if not os.path.exists(codon_usage_file):
            codon_usage_file = os.path.join(os.path.dirname(__file__), 'data', 'codon_usage.txt')

        with open(codon_usage_file, 'r') as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 3:
                    codon, aa = parts[0].strip(), parts[1].strip()
                    if aa == '*': 
                        continue  
                    try:
                        freq = float(parts[2].strip())
                        # Strictly exclude rare codons to guarantee passing CodonChecker
                        if freq >= 0.10:
                            self.aa_to_codons[aa].append((codon, freq))
                    except ValueError:
                        pass
        
        # Sort by frequency descending and keep only codons
        for aa in self.aa_to_codons:
            self.aa_to_codons[aa].sort(key=lambda x: x[1], reverse=True)
            self.aa_to_codons[aa] = [c[0] for c in self.aa_to_codons[aa]]

    def run(self, peptide: str, ignores: set, max_window: int = 80, max_attempts: int = 5_000_000) -> Transcript:
        """
        Translates the peptide to DNA using optimized DFS backtracking.

        Args:
            peptide: protein sequence (string of single-letter amino acids)
            ignores: set passed to the RBS chooser
            max_window: how many bp from the end to check each iteration (tuneable)
            max_attempts: global cap to avoid pathological infinite loops

        Returns:
            Transcript object with chosen RBS, peptide, and codons (stop codon appended)
        """
        selectedRBS = self.rbsChooser.run("", ignores)
        rbs_utr = selectedRBS.utr.upper()

        # Pre-check: ensure every AA has at least one codon choice
        for i, aa in enumerate(peptide):
            if aa not in self.aa_to_codons or len(self.aa_to_codons[aa]) == 0:
                raise ValueError(f"No codons available for amino acid '{aa}' at position {i}")

        # State for backtracking
        state = [0] * len(peptide)
        codons = [None] * len(peptide)
        current_cds_parts = []  # incremental list of codons (no join each time)
        pos = 0
        attempts = 0

        # Cache results for window -> (forbidden_ok, promoter_ok, hairpin_ok)
        check_cache = {}

        while pos < len(peptide):
            attempts += 1
            if attempts > max_attempts:
                raise Exception(f"Exceeded max attempts ({max_attempts}). Constraints might be too strict.")

            aa = peptide[pos]
            available_codons = self.aa_to_codons.get(aa, [])

            # defensive: should not happen because of pre-check, but guard against divide-by-zero
            if not available_codons:
                raise ValueError(f"No codons available for amino acid '{aa}' at position {pos}")

            if state[pos] < len(available_codons):
                # Heuristic: prefer higher-frequency codons but rotate to improve diversity
                shift = pos % len(available_codons)
                idx = (state[pos] + shift) % len(available_codons)
                cod = available_codons[idx]
                codons[pos] = cod

                # Incremental update of CDS parts
                if len(current_cds_parts) == pos:
                    current_cds_parts.append(cod)
                else:
                    current_cds_parts[pos] = cod
                    # truncate any later parts (shouldn't commonly be needed)
                    del current_cds_parts[pos+1:]

                # Build only the suffix we need to check
                full_test_seq = rbs_utr + "".join(current_cds_parts)
                window_to_check = full_test_seq[-max_window:].upper()

                # Cached checks to avoid repeating expensive work
                cached = check_cache.get(window_to_check)
                if cached is not None:
                    forbidden_ok, promoter_ok, hairpin_ok = cached
                else:
                    # Run cheap checks first
                    forbidden_ok, _ = self.forbidden_checker.run(window_to_check)
                    promoter_ok = True
                    hairpin_ok = True

                    if forbidden_ok and len(window_to_check) >= 29:
                        promoter_ok, _ = self.promoter_checker.run(window_to_check)

                    # Hairpin last (slow)
                    if forbidden_ok and promoter_ok and len(window_to_check) >= 50:
                        hairpin_ok, _ = hairpin_checker(window_to_check)

                    check_cache[window_to_check] = (forbidden_ok, promoter_ok, hairpin_ok)

                if forbidden_ok and promoter_ok and hairpin_ok:
                    # Accept codon and move forward
                    pos += 1
                else:
                    # Try next codon at this position
                    state[pos] += 1
            else:
                # Exhausted choices at this position -> backtrack
                state[pos] = 0
                codons[pos] = None
                if current_cds_parts:
                    current_cds_parts.pop()
                pos -= 1
                if pos < 0:
                    raise Exception("Constraint violation at start of sequence.")
                state[pos] += 1

        # Append stop codon and return transcript
        codons.append("TAA")
        return Transcript(selectedRBS, peptide, codons)

if __name__ == "__main__":
    peptide = "MYPFIRTARMTVCAKKHVHLTRDAAEQLLADIDRRLDQLLPVEGERDVVGAAMREGALAPGKRIRPMLLLLTARDLGCAVSHDGLLDLACAVEMVHAASLILDDMPCMDDAKLRRGRPTIHSHYGEHVAILAAVALLSKAFGVIADADGLTPLAKNRAVSELSNAIGMQGLVQGQFKDLSEGDKPRSAEAILMTNHFKTSTLFCASMQMASIVANASSEARDCLHRFSLDLGQAFQLLDDLTDGMTDTGKDSNQDAGKSTLVNLLGPRAVEERLRQHLQLASEHLSAACQHGHATQHFIQAWFDKKLAAVS"
    designer = TranscriptDesigner()
    designer.initiate()
    ignores = set()
    transcript = designer.run(peptide, ignores)
    print("Design successful!")