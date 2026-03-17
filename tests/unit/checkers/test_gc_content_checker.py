import pytest
from genedesign.checkers.gc_content_checker import GCContentChecker

@pytest.fixture
def gc_checker():
    """Fixture to initialize the GC checker before each test."""
    checker = GCContentChecker()
    checker.initiate()
    return checker

def test_passing_sequence(gc_checker):
    # 50% GC content
    seq = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT"
    result, bad_region = gc_checker.run(seq)
    
    assert result is True
    assert bad_region is None

def test_failing_high_gc_global(gc_checker):
    # 100% GC content
    seq = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"
    result, bad_region = gc_checker.run(seq)
    
    assert result is False
    assert "Global GC out of bounds" in bad_region

def test_failing_low_gc_window(gc_checker):
    # 50 As and 50 Cs. 
    # Global GC is exactly 50% (passes), but local 50bp windows are 0% and 100% (fails).
    seq = ("A" * 50) + ("C" * 50)
    result, bad_region = gc_checker.run(seq)
    
    assert result is False
    assert "Local GC out of bounds" in bad_region
