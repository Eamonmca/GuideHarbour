from probeharbour.blast_backend import BlastHit, _parse_hit


def test_parse_blast_hit_retains_raw_alignment_evidence():
    hit = _parse_hit("ENSMUST0001\t100.0\t32\t1\t32\t5\t36\t1e-10\t64.0")

    assert hit == BlastHit("ENSMUST0001", 100.0, 32, 1, 32, 5, 36, 1e-10, 64.0)
