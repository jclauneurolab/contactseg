rule download_template:
    params:
        template=config["template_flow"],
    output:
        template_txt = "resources/templateflow_template.txt"
    conda:
        "../envs/templateflow.yaml"
    script:
        "../scripts/template_flow.py"