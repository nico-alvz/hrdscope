from hrdscope.scars import Arm, Segment, hrd_loh, lst, merge_adjacent, ntai

ARMS = {"chr1": Arm(length=248_956_422, cen_start=121_700_000, cen_end=125_100_000)}
M = 1_000_000


def seg(start, end, major, minor, chrom="chr1"):
    return Segment(chrom, start * M, end * M, major, minor)


def test_merge_adjacent_joins_same_state():
    segs = [seg(0, 10, 1, 1), seg(10, 20, 1, 1), seg(20, 30, 2, 0)]
    assert [s.length for s in merge_adjacent(segs)] == [20 * M, 10 * M]


def test_loh_counts_large_interstitial_regions_only():
    segs = [seg(0, 20, 2, 0), seg(20, 30, 1, 1), seg(30, 40, 2, 0), seg(40, 248, 1, 1)]
    assert hrd_loh(segs) == 1  # 20 Mb region counts, 10 Mb region does not


def test_loh_excludes_whole_chromosome():
    assert hrd_loh([seg(0, 120, 2, 0), seg(120, 248, 3, 0)]) == 0


def test_ntai_counts_telomeric_imbalance_not_crossing_centromere():
    segs = [seg(0, 50, 2, 1), seg(50, 200, 1, 1), seg(200, 248, 3, 0)]
    assert ntai(segs, ARMS) == 2


def test_ntai_ignores_imbalance_crossing_centromere():
    segs = [seg(0, 130, 2, 1), seg(130, 248, 1, 1)]
    assert ntai(segs, ARMS) == 0


def test_lst_counts_breaks_between_large_segments_per_arm():
    # p arm: 0-121.7: two 50 Mb segments -> 1 break; q arm: 125.1-248: three segments -> 2 breaks
    segs = [seg(0, 50, 1, 1), seg(50, 121, 2, 1), seg(121, 126, 2, 1),
            seg(126, 160, 1, 1), seg(160, 200, 2, 2), seg(200, 248, 1, 0)]
    assert lst(segs, ARMS) == 3


def test_lst_smooths_out_small_segments():
    segs = [seg(0, 60, 1, 1), seg(60, 61, 3, 1), seg(61, 121, 1, 1)]
    assert lst(segs, ARMS) == 0
