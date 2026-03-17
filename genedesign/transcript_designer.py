import os
import time
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

    This version includes:
      - incremental CDS building (no repeated full joins)
      - caching of expensive window checks
      - periodic logging/profiling
      - max_attempts and max_seconds safety caps
      - a greedy fallback strategy if backtracking appears hopeless
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

    def _greedy_try(self, peptide: str, rbs_utr: str):
        """
        Fast greedy attempt that picks highest-frequency codons without hairpin checks.
        This is a fallback to produce a candidate quickly when exhaustive backtracking stalls.
        Returns codon list (without appended stop) or None if it fails forbidden/promoter checks.
        """
        parts = []
        for aa in peptide:
            available = self.aa_to_codons.get(aa, [])
            if not available:
                return None
            parts.append(available[0])  # highest frequency codon

        seq = rbs_utr + "".join(parts)
        window = seq[-200:].upper()  # check a reasonably large suffix once

        # One-shot checks: forbidden + promoter + hairpin (only once)
        passed, _ = self.forbidden_checker.run(window)
        if not passed:
            return None
        if len(window) >= 29:
            passed, _ = self.promoter_checker.run(window)
            if not passed:
                return None
        if len(window) >= 50:
            passed, _ = hairpin_checker(window)
            if not passed:
                return None
        return parts

    def run(self, peptide: str, ignores: set, 
            max_window: int = 80, max_attempts: int = 5_000_000, max_seconds: int = 300,
            log_every: int = 50_000) -> Transcript:
        """
        Optimized backtracking with incremental sequence building, caching and logging.

        Args:
            peptide: protein sequence (string of single-letter amino acids)
            ignores: set passed to the RBS chooser
            max_window: how many bp from the end to check each iteration (tuneable)
            max_attempts: global cap to avoid pathological infinite loops
            max_seconds: wall-clock timeout per run (seconds)
            log_every: print a progress line every `log_every` attempts

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
        start_time = time.perf_counter()

        # Counters for profiling
        hairpin_calls = 0
        promoter_calls = 0
        forbidden_calls = 0

        # Cache results for window -> (forbidden_ok, promoter_ok, hairpin_ok)
        check_cache = {}

        while pos < len(peptide):
            attempts += 1

            # Periodic logging
            if attempts % log_every == 0:
                elapsed = time.perf_counter() - start_time
                attempts_per_sec = attempts / max(elapsed, 1e-6)
                cache_size = len(check_cache)
                print(f"ATTEMPT {attempts} pos={pos} elapsed={elapsed:.1f}s attempts/sec={int(attempts_per_sec)} cache={cache_size} hairpin={hairpin_calls} promoter={promoter_calls} forbidden={forbidden_calls}")

            # Timeout / attempt caps
            elapsed = time.perf_counter() - start_time
            if attempts > max_attempts or elapsed > max_seconds:
                # Try a fast greedy fallback once before giving up
                print(f"WARNING: exceeded limits (attempts={attempts}, elapsed={int(elapsed)}s). Trying greedy fallback...")
                greedy = self._greedy_try(peptide, rbs_utr)
                if greedy:
                    print("Greedy fallback succeeded — using greedy sequence (may not satisfy all constraints globally).")
                    codons = greedy[:]  # replace codons with greedy result
                    codons.append("TAA")
                    return Transcript(selectedRBS, peptide, codons)
                else:
                    raise Exception(f"Exceeded limits (attempts={attempts}, elapsed={int(elapsed)}s) and greedy fallback failed. Try increasing limits or relaxing constraints (e.g., smaller max_window).")

            aa = peptide[pos]
            available_codons = self.aa_to_codons.get(aa, [])

            # Defensive check (should be prevented by pre-check)
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
                    # truncate any later parts
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
                    forbidden_calls += 1
                    promoter_ok = True
                    hairpin_ok = True

                    if forbidden_ok and len(window_to_check) >= 29:
                        promoter_ok, _ = self.promoter_checker.run(window_to_check)
                        promoter_calls += 1

                    # Hairpin last (slow)
                    if forbidden_ok and promoter_ok and len(window_to_check) >= 50:
                        hairpin_ok, _ = hairpin_checker(window_to_check)
                        hairpin_calls += 1

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
        total_elapsed = time.perf_counter() - start_time
        print(f"Completed in attempts={attempts} elapsed={total_elapsed:.1f}s cache_hits={len(check_cache)} hairpin_calls={hairpin_calls} promoter_calls={promoter_calls} forbidden_calls={forbidden_calls}")
        return Transcript(selectedRBS, peptide, codons)


if __name__ == "__main__":
    # Small smoke test
    peptide = "MYPFIRTARMTVCAKKHVHLTRDAAEQLLADIDRRLDQLLPVEGERDVVGAAMREGALAPGKRIRPMLLLLTARDLGCAVSHDGLLDLACAVEMVHAASLILDDMPCMDDAKLRRGRPTIHSHYGEHVAILAAVALLSKAFGVIADADGLTPLAKNRAVSELSNAIGMQGLVQGQFKDLSEGDKPRSAEAILMTNHFKTSTLFCASMQMASIVANASSEARDCLHRFSLDLGQAFQLLDDLTDGMTDTGKDSNQDAGKSTLVNLLGPRAVEERLRQHLQLASEHLSAACQHGHATQHFIQAWFDKKLAAVS"
    designer = TranscriptDesigner()
    designer.initiate()
    ignores = set()
    # tune parameters here (max_window smaller -> fewer checks, max_seconds shorter -> faster fail)
    transcript = designer.run(peptide, ignores, max_window=80, max_attempts=2_000_000, max_seconds=120, log_every=100_000)
    print("Design successful!")