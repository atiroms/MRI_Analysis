import os
def create_key(template, outtype=('nii.gz',), annotation_classes=None):
    if template is None or not template:
        raise ValueError('Template must be a valid format string')
    return template, outtype, annotation_classes
def infotodict(seqinfo):
    """Heuristic evaluator for determining which runs belong where
    allowed template fields - follow python string module:
    item: index within category
    subject: participant id
    seqitem: run number during scanning
    subindex: sub index within group
    """
    t1w = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_T1w')

    func_rest = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_bold')

    fieldmap = create_key('sub-{subject}/{session}/fmap/sub-{subject}_{session}_fieldmap')
    info = {t1w: [], func_rest: [], fieldmap: []}
    #info = {t1w: [], fieldmap: []}

    for idx, s in enumerate(seqinfo):
        if s.dim1==256 and s.dim2==256 and s.dim3==200:
            info[t1w].append(s.series_id)
        if s.dim1==64 and s.dim2==64 and s.dim3==10000:
            info[func_rest].append(s.series_id)
        if s.dim1==64 and s.dim2==64 and s.dim3==160 and (s.series_description=='+Fieldmap_SBPRS' or s.series_description=='+Fieldmap'):
            info[fieldmap].append(s.series_id)
    return info
