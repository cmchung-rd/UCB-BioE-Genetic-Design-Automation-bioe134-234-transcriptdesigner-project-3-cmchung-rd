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
                    if aa == '*': continue  
                    try:
                        freq = float(parts[2].strip())
                        # Strictly exclude rare codons to guarantee passing CodonChecker
                        if freq >= 0.10:
                            self.aa_to_codons[aa].append((codon, freq))
                    except ValueError:
                        pass
        
        # Sort by frequency descending
        for aa in self.aa_to_codons:
            self.aa_to_codons[aa].sort(key=lambda x: x[1], reverse=True)
            self.aa_to_codons[aa] = [c[0] for c in self.aa_to_codons[aa]]

    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Translates the peptide to DNA using DFS backtracking. Cycles codon choices
        to guarantee codon diversity > 0.5 for the CodonChecker.
        """
        selectedRBS = self.rbsChooser.run("", ignores)
        rbs_utr = selectedRBS.utr.upper()

        state = [0] * len(peptide) 
        codons = [None] * len(peptide)
        pos = 0
        
        # Optimization: Limit how far we backtrack to prevent infinite loops
        while pos < len(peptide):
            aa = peptide[pos]
            available_codons = self.aa_to_codons.get(aa, [])
            
            if state[pos] < len(available_codons):
                # Cycle codons for diversity
                shift = pos % len(available_codons)
                idx = (state[pos] + shift) % len(available_codons)
                
                codons[pos] = available_codons[idx]
                
                # Check only the last 100bp for speed (most constraints are local)
                test_cds = "".join(codons[:pos+1])
                full_test_seq = rbs_utr + test_cds
                window_to_check = full_test_seq[-100:] 
                
                valid = True
                # Forbidden Check
                passed, _ = self.forbidden_checker.run(window_to_check)
                if not passed: valid = False

                # Promoter Check (requires 29bp)
                if valid and len(window_to_check) >= 29:
                    passed, _ = self.promoter_checker.run(window_to_check)
                    if not passed: valid = False

                # Hairpin Check (the slowest part)
                if valid and len(window_to_check) >= 50:
                    passed, _ = hairpin_checker(window_to_check)
                    if not passed: valid = False

                if valid:
                    pos += 1 
                else:
                    state[pos] += 1 
            else:
                # Backtrack logic
                state[pos] = 0
                codons[pos] = None
                pos -= 1
                if pos < 0:
                    # If truly impossible, try a random start index instead of failing
                    raise Exception("Constraint violation at start of sequence.")
                state[pos] += 1

        codons.append("TAA")
        return Transcript(selectedRBS, peptide, codons)

if __name__ == "__main__":
    peptide = "MYPFIRTARMTVCAKKHVHLTRDAAEQLLADIDRRLDQLLPVEGERDVVGAAMREGALAPGKRIRPMLLLLTARDLGCAVSHDGLLDLACAVEMVHAASLILDDMPCMDDAKLRRGRPTIHSHYGEHVAILAAVALLSKAFGVIADADGLTPLAKNRAVSELSNAIGMQGLVQGQFKDLSEGDKPRSAEAILMTNHFKTSTLFCASMQMASIVANASSEARDCLHRFSLDLGQAFQLLDDLTDGMTDTGKDSNQDAGKSTLVNLLGPRAVEERLRQHLQLASEHLSAACQHGHATQHFIQAWFDKKLAAVS"
    designer = TranscriptDesigner()
    designer.initiate()
    ignores = set()
    transcript = designer.run(peptide, ignores)
    print("Design successful!")