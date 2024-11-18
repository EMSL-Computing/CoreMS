from corems.molecular_id.search.compoundSearch import LowResMassSpectralMatch


def get_rt_ri_pairs(gcms_ref_obj, sql_obj=None):
    lowResSearch = LowResMassSpectralMatch(
        gcms_ref_obj, sql_obj=sql_obj, calibration=True
    )

    lowResSearch.run()

    dict_ri_rt = {}

    list_of_compound_obj = {}

    for gcms_peak in gcms_ref_obj:
        # has a compound matched
        if gcms_peak:
            compound_obj = gcms_peak.highest_ss_compound

            if not compound_obj.ri in dict_ri_rt.keys():
                dict_ri_rt[compound_obj.ri] = [
                    (gcms_peak.mass_spectrum.retention_time, compound_obj)
                ]

            else:
                dict_ri_rt[compound_obj.ri].append(
                    (gcms_peak.mass_spectrum.retention_time, compound_obj)
                )
            if gcms_ref_obj.parameters.gc_ms.verbose_processing:
                print(
                    compound_obj.name,
                    gcms_peak.mass_spectrum.retention_time,
                    compound_obj.spectral_similarity_score,
                )

    ris = [i for i in dict_ri_rt.keys()]
    rts = [
        max(i, key=lambda c: c[1].spectral_similarity_score)[0]
        for i in dict_ri_rt.values()
    ]

    rt_ri_pairs = list(zip(rts, ris))

    return rt_ri_pairs
