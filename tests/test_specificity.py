from probeharbour.specificity import SpecificityEvidence, SpecificityRequest, SpecificityService


class FakePanelApi:
    def __init__(self):
        self.calls = 0

    def evaluate(self, request):
        self.calls += 1
        return SpecificityEvidence(request, target_raw_score=64, best_offtarget_raw_score=16, off_target_count=1, best_offtarget_id="off-target")


def test_panel_api_backend_is_cache_first_and_keeps_method_provenance():
    request = SpecificityRequest("ACTG", "tx1", "refseq", "2026-07", "blast-api-v2")
    backend = FakePanelApi()
    service = SpecificityService(backend)
    first, cached_first = service.evaluate(request)
    second, cached_second = service.evaluate(request)
    assert first.specificity_percent == 75.0
    assert second == first
    assert not cached_first and cached_second
    assert backend.calls == 1
    assert first.request.method_id == "blast-api-v2"
