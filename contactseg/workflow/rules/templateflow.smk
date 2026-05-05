rule download_template:
    params:
        template=config["template_flow"],
        atlas=config["template_atlas"],
        atlas_desc=config["atlas_desc"],
    input:
        default_template_path=str(
            Path(workflow.basedir).parent.parent
            / "resources/atlases/tpl-MNI152NLin2009cSym_res-1_T1w.nii.gz"
        ),
        default_atlas_path=str(
            Path(workflow.basedir).parent.parent
            / "resources/atlases/tpl-MNI152NLin2009cSym_res-1_atlas-CerebrA_dseg.nii"
        ),
        default_tsv_path=str(
            Path(workflow.basedir).parent.parent
            / "resources/atlases/tpl-MNI152NLin2009cSym_atlas-CerebA_dseg.tsv"
        ),
    output:
        template_txt="resources/templateflow_template.txt",
    conda:
        "../envs/templateflow.yaml"
    script:
        "../scripts/template_flow.py"
