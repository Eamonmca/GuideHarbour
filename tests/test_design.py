from probeharbour.design import generate_junction_candidates, select_candidates
from probeharbour.models import Provenance

def test_generate_all_valid_slides_and_keep_slide_evidence():
    probes = generate_junction_candidates("A" * 80, junction_position=40, arm_length=16, slides=[-2, 0, 2], provenance=Provenance("test", "fixture"), candidate_prefix="tx")
    assert [probe.shift_nt for probe in probes] == [-2, 0, 2]
    assert probes[1].on_junction

def test_selection_keeps_reasons_for_losing_candidates():
    decisions = select_candidates([{"candidate_id":"a","recalculated_score":1,"specificity_percent":90,"audit_ok":True},{"candidate_id":"b","recalculated_score":.9,"specificity_percent":100,"audit_ok":True}])
    assert decisions[0].selected and decisions[1].reasons == ("lower ranked eligible candidate",)

def test_selection_accepts_missing_specificity_before_blast():
    decisions = select_candidates([
        {"candidate_id": "lower", "recalculated_score": 0.5, "specificity_percent": None},
        {"candidate_id": "higher", "recalculated_score": 0.9, "specificity_percent": None},
    ])
    assert decisions[0].candidate_id == "higher"
    assert decisions[0].selected
