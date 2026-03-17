class GCContentChecker:
    """
    Checks the GC content of a given DNA sequence to ensure it falls 
    within acceptable synthesis thresholds (both globally and locally).
    """

    def __init__(self):
        self.min_gc = 0.0
        self.max_gc = 0.0
        self.window_size = 0

    def initiate(self) -> None:
        """Sets the acceptable bounds for GC content and the window size."""
        self.min_gc = 0.30  # 30% minimum GC content
        self.max_gc = 0.70  # 70% maximum GC content
        self.window_size = 50 # 50 bp sliding window

    def run(self, dna: str) -> tuple[bool, str]:
        """
        Checks if the DNA sequence meets GC content requirements.
        
        Parameters:
            dna (str): The DNA sequence to analyze.
            
        Returns:
            tuple: (True, None) if the sequence passes.
                   (False, problematic_string) if it fails.
        """
        if not dna:
            return False, "Empty sequence"

        dna = dna.upper()
        
        # 1. Check global GC content
        gc_count = dna.count('G') + dna.count('C')
        global_gc = gc_count / len(dna)
        if global_gc < self.min_gc or global_gc > self.max_gc:
            return False, f"Global GC out of bounds: {global_gc:.2f}"

        # 2. Check local GC content using a sliding window
        if len(dna) >= self.window_size:
            for i in range(len(dna) - self.window_size + 1):
                window = dna[i:i + self.window_size]
                window_gc = (window.count('G') + window.count('C')) / self.window_size
                
                if window_gc < self.min_gc or window_gc > self.max_gc:
                    return False, f"Local GC out of bounds in window: {window}"
                
        return True, None

if __name__ == "__main__":
    checker = GCContentChecker()
    checker.initiate()
    
    good_seq = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT"
    bad_seq =  "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"
    
    print(f"Good Seq Result: {checker.run(good_seq)}")
    print(f"Bad Seq Result: {checker.run(bad_seq)}")