import pytest

from probeharbour.coordinates import coordinate_map_from_exons


@pytest.mark.parametrize("strand, expected", [
    ("+", tuple(range(100, 110)) + tuple(range(200, 210))),
    ("-", tuple(range(209, 199, -1)) + tuple(range(109, 99, -1))),
])
def test_coordinate_mapping_respects_transcript_orientation(strand, expected):
    mapping = coordinate_map_from_exons("tx1", "chr1", strand, [(100, 110), (200, 210)])
    assert tuple(mapping.spliced_to_genomic(position) for position in range(mapping.length)) == expected
    assert all(mapping.genomic_to_spliced(position) == index for index, position in enumerate(expected))
    assert mapping.junction_positions == (10,)
    assert mapping.junction_id(10) == "tx1:splice:10"
