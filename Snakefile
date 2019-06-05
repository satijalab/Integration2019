seurat_code = "software/seurat"

rule preprocessing:
    input:
        dl_script = "preprocessing_scripts/{dataset}.sh",
        r_script = "preprocessing_scripts/{dataset}.R",
        seurat = "software/seurat"
    params:
        source = "orig" # switch to "dropbox" to download archived version from dropbox
    conda: 
        "envs/seurat.yaml"
    output:
        "seurat_objects/{dataset}.rds"
    shell:
        """
        sh {input.dl_script} {params.source}
        Rscript {input.r_script} {seurat_code} {output}
        """

rule no_integration:
    input:
        data = "seurat_objects/integrated_{dataset}.rds",
        script = "analysis_code/integration/nointegration.R"
    conda:
        "envs/seurat.yaml"
    output:
        "seurat_objects/all_{dataset}_no_integration.rds"
    shell:
        "Rscript {input.script} {seurat_code} {input.data} {output}"

rule no_integration_all:
    input:
        expand("seurat_objects/all_{dataset}_no_integration.rds", dataset = ["pancreas", "bipolar"])

rule scrub_ica:
    input:
        script = "analysis_code/citeseq/scrub_ica.py",
        data = "seurat_objects/ica_bone_marrow.rds"		
    conda:
        "envs/scrublet.yaml"
    output:
        expand("raw_data/ica/{dataset}_doublets.tsv", dataset = ["MantonBM1", "MantonBM2", "MantonBM3", "MantonBM4", "MantonBM5", "MantonBM6", "MantonBM7", "MantonBM8"])
    shell:
        "python3 {input.script} --path raw_data/ica/"

###################### Full Integrations ###########################

rule integrate_pancreas:
    input:
        data = expand("seurat_objects/{dataset}.rds", dataset = ["celseq", "celseq2", "fluidigmc1", "smartseq2", "indrop"]),
        script = "analysis_code/integration/integrate_pancreas.R"
    conda: 
        "envs/seurat.yaml"
    output:
        "seurat_objects/integrated_pancreas.rds",
        "tables/Supplementary_Table_2.csv"
    shell:
        "Rscript {input.script} {seurat_code} {input.data} {output}"

rule integrate_bipolar:
    input:
        data = expand("seurat_objects/{dataset}.rds", dataset = ["bipolar"]),
        script = "analysis_code/integration/integrate_bipolar.R"
    conda:
        "envs/seurat.yaml"
    output:
        "seurat_objects/integrated_bipolar.rds"
    shell:
        "Rscript {input.script} {input.data} {output} {seurat_code}"

rule integrate_mca:
    input:
        data = expand("seurat_objects/{dataset}.rds", dataset = ["mca_tm_droplet", "mca_tm_facs"]),
        script = "analysis_code/integration/integrate_mca.R"
    output:
        "seurat_objects/integrated_mca.rds"
    shell:
        "Rscript {input.script} {input.data} {output}"

rule integrate_hca_bone:
    input:
        data = "seurat_objects/ica_bone_marrow.rds",
        script = "analysis_code/integration/integrate_bone_marrow.R",
        scrublet_results = expand("raw_data/ica/ICA_scrublet_results/{dataset}_doublets.tsv", dataset = ["MantonBM1", "MantonBM2", "MantonBM3", "MantonBM4", "MantonBM5", "MantonBM6", "MantonBM7", "MantonBM8"])
    output:
        "seurat_objects/integrated_bone_marrow.rds"
    shell:
        "Rscript {input.script} {input.data} {output}"

########################## Holdout Integrations  ##############################

rule generate_integration_holdouts:
    input:
        data = "seurat_objects/integrated_{dataset}.rds",
        script = "analysis_code/integration/{dataset}_celltype_holdouts.R"
    conda:
        "envs/seurat.yaml"
    output:
        "tables/integrated_{dataset}_celltype_holdout_table.csv",
        "seurat_objects/integrated_{dataset}_celltype_holdouts_none.rds"
    shell:
        "Rscript {input.script} {input.data} {output} {seurat_code}"

rule integrate_holdouts_seuratV3:
    input:
        data = "seurat_objects/integrated_{dataset}_celltype_holdouts_none.rds",
        script = "analysis_code/integration/integrate_{dataset}_celltype_holdouts_seuratV3.R"
    conda:
        "envs/seurat.yaml"
    output:
        "seurat_objects/integrated_{dataset}_celltype_holdouts_seuratV3.rds"
    shell:
        "Rscript {input.script} {input.data} {output} {seurat_code}"

rule integrate_holdouts_seuratV2:
    input:
        data = "seurat_objects/integrated_{dataset}_celltype_holdouts_none.rds",
        script = "analysis_code/integration/integrate_{dataset}_celltype_holdouts_seuratV2.R",
    conda:
        "envs/seurat2.yaml"
    output:
        "seurat_objects/integrated_{dataset}_celltype_holdouts_seuratV2.rds"
    shell:
        """
        Rscript {input.script} {input.data} {output} {seurat_code}
        """

rule integrate_holdouts_mnncorrect:
    input:
        data = "seurat_objects/integrated_{dataset}_celltype_holdouts_none.rds",
        script = "analysis_code/integration/integrate_{dataset}_celltype_holdouts_mnnCorrect.R"
    conda:
        "envs/scran.yaml"
    output:
        "seurat_objects/integrated_{dataset}_celltype_holdouts_mnnCorrect.rds"
    shell:
        """
        Rscript {input.script} {input.data} {output} {seurat_code}
        """

rule integrate_holdouts_scanorama:
    input:
        data = "seurat_objects/integrated_{dataset}_celltype_holdouts_none.rds",
        write_txt = "analysis_code/integration/{dataset}_scanorama_write_txt.R",
        convert_txt = "analysis_code/integration/scanorama_process.py",
        py_script = "analysis_code/integration/integrate_celltype_holdouts_scanorama.py",
        r_script = "analysis_code/integration/integrate_{dataset}_celltype_holdouts_scanorama.R",
        software = "software/scanorama/"
    conda:
        "envs/scanorama.yaml",
    output:
        "seurat_objects/integrated_{dataset}_celltype_holdouts_scanorama.rds"
    shell:
        """
        Rscript {input.write_txt} {input.data} {seurat_code}
        python {input.convert_txt} analysis_data/{wildcards.dataset}/{wildcards.dataset}_conf.txt
        python {input.py_script} analysis_data/{wildcards.dataset}/{wildcards.dataset}_conf2.txt
        Rscript {input.r_script} {input.data} {output} {seurat_code}
        """

rule integrate_holdouts_scmerge:
    input:
        data = "seurat_objects/integrated_{dataset}_celltype_holdouts_none.rds",
        software = "software/scMerge/",
        script = "analysis_code/integration/integrate_{dataset}_celltype_holdouts_scmerge.R"
    output:
        "seurat_objects/integrated_{dataset}_celltype_holdouts_scmerge.rds"
    shell:
        "Rscript {input.script} {input.software} {input.data} {output}"

rule integrate_holdouts_all:
    input:
        expand("seurat_objects/integrated_{dataset}_celltype_holdouts_seuratV3.rds", dataset = ["pancreas", "bipolar"]),
        expand("seurat_objects/integrated_{dataset}_celltype_holdouts_seuratV2.rds", dataset = ["pancreas", "bipolar"]),
        expand("seurat_objects/integrated_{dataset}_celltype_holdouts_mnnCorrect.rds", dataset = ["pancreas", "bipolar"]),
        expand("seurat_objects/integrated_{dataset}_celltype_holdouts_scanorama.rds", dataset = ["pancreas", "bipolar"])


########################## Figures ###################################

rule integration_dimplots:
    input:
        "figure_code/integration_dimplots.R",
        expand("seurat_objects/integrated_pancreas_celltype_holdouts_{method}.rds", method = ["seuratV3", "seuratV2", "mnnCorrect", "scanorama"])
    conda:
        "envs/figures.yaml"
    output:
        "figures/integration_dimplots.rda"
    shell:
        "Rscript {input} {output} {seurat_code}"

rule anchor_score_figure:
    input:
        data = expand("seurat_objects/{dataset}.rds", dataset = ["celseq", "celseq2", "smartseq2", "fluidigmc1", "indrop", "integrated_pancreas"]),
        script = "figure_code/anchor_scores.R"
    conda:
        "envs/figures.yaml"
    output:
        "figures/anchor_scores.pdf",
        "figures/anchor_scores.rds",
        "figures/anchor_barplot.pdf",
        "figures/anchor_barplot.rds"
    shell:
        "Rscript {input.script} {input.data} {output} {seurat_code}"

rule metric_figures:
    input:
        data = expand("analysis_data/integrated_{dataset}_celltype_holdouts_{method}_metrics.rds",
            method = ["seuratV3", "seuratV2", "mnnCorrect", "scanorama", "none"], dataset = ["pancreas", "bipolar"]),
        script = "figure_code/metric_figures.R"
    conda:
        "envs/figures.yaml"
    output:
        "figures/silhouette.rds",
        "figures/mixing_ls_metrics.rds",
    shell:
        "Rscript {input.script} {input.data} {output} {seurat_code}"

rule figure2:
    input:
        "figure_code/figure2.R",
        "figures/integration_dimplots.rda",
        "figures/anchor_scores.rds",
        "figures/anchor_barplot.rds",
        "figures/silhouette.rds",
        "figures/mixing_ls_metrics.rds",
    conda:
        "envs/figures.yaml"
    output:
        "figures/supplementary/pancreas_full_umaps.pdf",
        "figures/figure2.pdf"
    shell:
        "Rscript {input} {output} {seurat_code}"


rule figure3_holdouts:
    input:
        "figure_code/figure3_holdouts.R",
        "analysis_data/pancreas/indrop4-alpha-projection-results.rds",
        "analysis_data/pancreas/all_accuracy.rds",
        "analysis_data/bipolar/all_accuracy.rds"
    conda:
        "envs/figures.yaml"
    output:
        "figures/figure3_confusion_matrices.pdf",
        "figures/figure3_accuracy_boxplot.pdf",
        "figures/figure3_score_histograms.pdf",
        "tables/Supplementary_Table_3.csv"
    shell:
        "Rscript {input} {output} {seurat_code}"

rule figure3_coembed_umaps:
    input:
        "figure_code/pfc_atac_coembed_umaps.R",
        "seurat_objects/pfc_atac_coembed.rds"
    output:
        "figures/figure4_pfc_atac_coembed_umaps.pdf"
    shell:
        "Rscript {input} {output}"

rule citeseq_cv_fig:
    input:
        "figure_code/citeseq_cv.R",
        "seurat_objects/citeseq.rds",
        "seurat_objects/citeseq_crossvalidation.rds"
    output:
        "figures/citeseq/cv_correlation.png",
        "figures/citeseq/cv_correlation_select.png",
        "figures/citeseq/cv_correlation_select.rds"
    shell:
        "Rscript {input} {output}"


rule citeseq_ds_fig:
    input:
        "figure_code/citeseq_downsampling_figures.R",
        "seurat_objects/citeseq.rds",
        "seurat_objects/citeseq_downsampling/1000_feature_downsampling.rds"
    output:
        "figures/citeseq/nfeature_correlation.pdf",
        "figures/citeseq/nfeature_select_correlation.pdf",
        "figures/citeseq/nfeature_heatmap_correlation.pdf",
        "figures/citeseq/nfeature_heatmap_correlation.rds"
    shell:
        "Rscript {input} {output}"

rule mca_figure:
    input:
        "figure_code/mca_figure.R",
        "seurat_objects/integrated_mca.rds"
    output:
        "figures/figureS3.pdf"
    shell:
        "Rscript {input} {output}"

rule bipolar_holdout_tsnes:
    input:
        "figure_code/supplementary/bipolar_holdout_tsnes.R",
        expand("seurat_objects/integrated_bipolar_celltype_holdouts_{method}.rds", method = ["seuratV3", "seuratV2", "mnnCorrect", "scanorama"])
    output:
        "figures/supplementary/bipolar_holdout_tsnes.pdf"
    shell:
        "Rscript {input} {output}"

rule all_figures:
    input:
        "figures/figure4.pdf",
        "figures/figure2.pdf"

########################## Metrics ###################################

rule integration_metrics:
    input:
        data = "seurat_objects/integrated_{dataset}_celltype_holdouts_{method}.rds",
        script = "analysis_code/integration/integration_metrics.R"
    conda:
        "envs/seurat.yaml"
    output:
        "analysis_data/integrated_{dataset}_celltype_holdouts_{method}_metrics.rds"
    shell:
        "Rscript {input.script} {input.data} {output} {seurat_code}"

########################## Projections  ##############################

rule pancreas_holdout_references:
    input:
        "seurat_objects/celseq.rds",
        "seurat_objects/celseq2.rds",
        "seurat_objects/smartseq2.rds",
        "seurat_objects/fluidigmc1.rds",
        "seurat_objects/indrop.rds",
        "seurat_objects/integrated_pancreas.rds"
    params:
        holdout="{holdout}",
        query="{query}"
    conda:
        "envs/seurat.yaml"
    output:
        "analysis_data/pancreas/reference-without-{query}-{holdout}.rds",
        "analysis_data/pancreas/query-without-{query}-{holdout}.rds"
    shell:
        "Rscript ./analysis_code/projection/pancreas_generate_holdout_references.R {input} {params.holdout} {params.query} {seurat_code}"

rule project_pancreas_holdouts:
    input:
        "analysis_code/projection/holdout_accuracy.R",
        "analysis_data/pancreas/reference-without-{query}-{holdout}.rds",
        "analysis_data/pancreas/query-without-{query}-{holdout}.rds",
    params:
        holdout="{holdout}",
        query="{query}"
    conda:
        "envs/projection.yaml"
    output:
        "analysis_data/pancreas/{query}-{holdout}-accuracy.rds",
        "analysis_data/pancreas/{query}-{holdout}-projection-results.rds"
    shell:
        "Rscript {input} {params.holdout} {params.query} {output} {seurat_code}"

rule pancreas_accuracy:
    input:
        expand("analysis_data/pancreas/{query}-{holdout}-accuracy.rds",
            holdout = ['alpha', 'acinar', 'activated_stellate', 'beta', 'delta', 'ductal', 'endothelial', 'gamma', 'quiescent_stellate'],
            query = ['celseq', 'celseq2', 'smartseq2', 'fluidigmc1', 'indrop1', 'indrop2', 'indrop3', 'indrop4'])
    conda:
        "envs/seurat.yaml"
    output:
        "analysis_data/pancreas/all_accuracy.rds"
    shell:
        "Rscript ./analysis_code/projection/collate_accuracy.R {output} {input}"

rule bipolar_holdout_references:
    input:
        "seurat_objects/bipolar.rds",
    params:
        holdout="{holdout}",
        query="{query}"
    conda:
        "envs/seurat.yaml"
    output:
        "analysis_data/bipolar/reference-without-{query}-{holdout}.rds",
        "analysis_data/bipolar/query-without-{query}-{holdout}.rds"
    shell:
        "Rscript ./analysis_code/projection/bipolar_generate_holdout_references.R {input} {params.holdout} {params.query} {seurat_code}"

rule project_bipolar_holdouts:
    input:
        "analysis_code/projection/holdout_accuracy.R",
        "analysis_data/bipolar/reference-without-{query}-{holdout}.rds",
        "analysis_data/bipolar/query-without-{query}-{holdout}.rds"
    params:
        holdout="{holdout}",
        query="{query}"
    output:
        "analysis_data/bipolar/{query}-{holdout}-accuracy.rds",
        "analysis_data/bipolar/{query}-{holdout}-projection-results.rds"
    shell:
        "Rscript {input} {params.holdout} {params.query} {output} {seurat_code}"

rule bipolar_accuracy:
    input:
        expand("analysis_data/bipolar/{query}-{holdout}-accuracy.rds",
            holdout = ['Amacrine_cells', 'Muller_glia', 'BC1A', 'BC1B', 'BC2', 'BC3A', 'BC3B', 'BC4', 'BC5A', 'BC5B', 'BC5C', 'BC5D', 'BC6', 'BC7', 'BC8_BC9','Cone_photoreceptors', 'RBC', 'Rod_photoreceptors' ],
            query=['1', '2', '3', '4', '5', '6'])
    conda:
        "envs/seurat.yaml"
    output:
        "analysis_data/bipolar/all_accuracy.rds"
    shell:
        "Rscript ./analysis_code/projection/collate_accuracy.R {output} {input}"

rule all_accuracy:
    input:
        "analysis_data/bipolar/all_accuracy.rds",
        "analysis_data/pancreas/all_accuracy.rds"

###################### Spatial  ###########################

rule preprocess_starmap:
    input:
        dl_script = "preprocessing_scripts/starmap.sh",
        python_script = "analysis_code/spatial/get_spatial_positions.py",
        r_script = "preprocessing_scripts/starmap.R" 
    params:
        source = "orig"
    conda:
        "envs/seurat.yaml"
    output:
        "seurat_objects/20180505_BY3_1kgenes.rds",
        "seurat_objects/20180410-BY3_1kgenes.rds"
    shell:
        """
        sh {input.dl_script} {params.source}
        python {input.python_script} --starmap_path=raw_data/spatial/starmap/ --data_path=raw_data/spatial/starmap/
        Rscript {input.r_script} {output} {seurat_code}
        """

rule get_imputation_files:
    input:
        expand("seurat_objects/{dataset}.rds",
        dataset = ["osm_fish", "allen_brain", "20180505_BY3_1kgenes", "dropseq_cortex"])


rule impute_starmap:
    input:
        "analysis_code/spatial/impute_starmap.R",
        "seurat_objects/20180505_BY3_1kgenes.rds",
        "seurat_objects/20180410-BY3_1kgenes.rds",
        "seurat_objects/allen_brain.rds",
        "seurat_objects/dropseq_cortex.rds"
    conda:
        "envs/seurat.yaml"
    output:
        "seurat_objects/20180505_BY3_1kgenes_imputed.rds",
        "seurat_objects/20180410-BY3_1kgenes_imputed.rds", 
        "seurat_objects/integrated_starmap.rds"
    shell:
        "Rscript {input} {output} {seurat_code}"

rule spatial_crossvalidation:
    input:
        "seurat_objects/20180505_BY3_1kgenes.rds",
        "seurat_objects/allen_brain.rds",
        "seurat_objects/dropseq_cortex.rds"
    output:
        "seurat_objects/spatial_leaveout/starmap_1_Cux2.rds",
        "analysis_data/starmap_markers.tsv"
    shell:
        "Rscript analysis_code/spatial/starmap_leaveout.R"

# fig6
rule analyze_imputed_spatial:
    input:
        "seurat_objects/osm_fish.rds",
        "seurat_objects/20180505_BY3_1kgenes_imputed.rds",
        "seurat_objects/integrated_starmap.rds",
        "seurat_objects/dropseq_cortex.rds",
        "seurat_objects/spatial_leaveout/starmap_1_Cux2.rds"
    output:
        "figures/spatial/starmap_rep1_vs_rep2_moran.png"
    shell:
        "Rscript analysis_code/spatial/analyze_imputed_starmap.R"

# supplementary note 2

rule sn2a:
    input:
        "analysis_code/spatial/starmap_downsampling_comparison.R",
        "seurat_objects/20180505_BY3_1kgenes.rds",
        "seurat_objects/20180410-BY3_1kgenes.rds",
        "seurat_objects/integrated_starmap.rds",
        "seurat_objects/allen_brain.rds",
        "seurat_objects/dropseq_cortex.rds",
        "analysis_data/starmap_markers.tsv"
    conda:
        "envs/seurat.yaml"
    output:
        "figures/sn2/starmap_downsample_guided_cluster_markers_ss2.png",
        "figures/sn2/starmap_downsample_ss2_random_vs_guided_boxplots.png",
        "figures/sn2/starmap_downsample_ss2_random_boxplots.png"
    shell:
        "Rscript {input} {output} {seurat_code}" 

rule sn2b:
    input:
        "analysis_code/spatial/redundancy_vs_accuracy.R",
        "seurat_objects/20180505_BY3_1kgenes.rds",
        "seurat_objects/20180410-BY3_1kgenes.rds",
        "seurat_objects/integrated_starmap.rds",
        "seurat_objects/allen_brain.rds",
        "seurat_objects/dropseq_cortex.rds"
    conda:
        "envs/seurat.yaml"
    output:
        "figures/sn2/starmap_downsample_guided_cluster_markers_ss2.png",
        "figures/sn2/starmap_downsample_ss2_random_vs_guided_boxplots.png",
        "figures/sn2/starmap_downsample_ss2_random_boxplots.png"
    shell:
        "Rscript {input} {output} {seurat_code}" 

###################### Protein  ###########################

rule make_citeseq:
    input:
        "seurat_objects/citeseq.rds"

rule make_hca:
    input:
        "seurat_objects/ica_bone_marrow.rds"

# 4A
rule crossvalidate_citeseq:
    input:
        "analysis_code/citeseq/cite_cross_validations.R",
        "seurat_objects/citeseq.rds"
    conda:
        "envs/seurat.yaml"
    output:
        "seurat_objects/citeseq_crossvalidation.rds"
    shell:
        "Rscript {input} {output} {seurat_code}"

# 4B
rule feature_downsample:
    input:
        "seurat_objects/citeseq.rds"
    conda:
        "envs/seurat.yaml"
    output:
        "seurat_objects/citeseq_downsampling/1000_feature_downsampling.rds"
    shell:
        "Rscript analysis_code/citeseq/cite_feature_downsampling.R {input} {seurat_code}"

rule impute_citeseq:
    input:
        "analysis_code/citeseq/impute_citeseq.R",
        "seurat_objects/citeseq.rds",
        "seurat_objects/integrated_bone_marrow.rds"
    conda:
        "envs/seurat.yaml"
    output:
        "seurat_objects/imputed_citeseq.rds"
    shell:
        "Rscript {input} {output} {seurat_code}"

# 4CDF
rule analyze_imputed_citeseq:
    input:
        "seurat_objects/imputed_citeseq.rds",
    output:
        "analysis_data/cd69_markers.tsv"
    shell:
        "Rscript analysis_code/citeseq/analyze_imputed_citeseq.R"

# supplement
rule go_enrichment:
    input:
        "analysis_data/cd69_markers.tsv"
    output:
        "figures/citeseq/go_terms_cd69_bp.pdf"
    shell:
        "Rscript analysis_code/citeseq/enrichment_analysis.R"

# E
rule download_facs_data:
    output:
        "raw_data/bone_marrow/CD69.umi.txt"
    shell:
        "sh preprocessing_scripts/cd69_facs.sh"

rule facs_analysis:
    input:
        "analysis_data/cd69_markers.tsv",
        "raw_data/bone_marrow/CD69.umi.txt"
    output:
        "figures/citeseq/facs_volcano.png"
    shell:
        "Rscript analysis_code/citeseq/analyze_facs.R"

###################### ATAC ###########################

rule coembed_atac:
    input:
        script="analysis_code/atac/coembed_pfc.R",
        rna="seurat_objects/allen_brain_upper.rds",
        atac="seurat_objects/atac_pfc.rds"
    output:
        coembed="seurat_objects/pfc_atac_coembed.rds",
        predictions="seurat_objects/pfc_atac_predictions.rds",
        pv="analysis_data/atac/clusters/pvalb.txt"
    shell:
        "Rscript {input.script} {input.rna} {input.atac} {output.coembed} {output.predictions}"

rule plot_atac_coembed:
    input:
        "seurat_objects/pfc_atac_coembed.rds"
    output:
        "figures/atac/pfc/coembed_atac_dataset.png",
        "figures/atac/pfc/coembed_atac_celltype.png",
        "figures/atac/pfc/coembed_atac_celltype_only_atac.png"
    shell:
        "Rscript figure_code/atac_coembed_umap.R"

rule find_da_peaks:
    input:
        "seurat_objects/pfc_atac_predictions.rds"
    output:
        "analysis_data/atac/da_peaks/pv.bed"
    shell:
        "Rscript analysis_code/atac/DA_peaks.R {input}"

rule find_motifs:
    input:
        "analysis_data/atac/da_peaks/pv.bed"
    output:
        "analysis_data/atac/homer_pv/homerResults.html"
    shell:
        "sh analysis_code/atac/find_motifs.sh"

rule split_clusters:
    input:
        "analysis_data/atac/clusters/pvalb.txt"
    output:
        "analysis_data/atac/clusters/pvalb.bw"
    conda:
        "envs/atac.yaml"
    shell:
        "sh analysis_code/atac/split_clusters.sh"


rule plot_coverage:
    input:
        "analysis_data/atac/clusters/pvalb.bw"
    output:
        "figures/atac/pfc/gad2.png"
    shell:
        "Rscript analysis_code/atac/plot_tracks.R"

### 10x ATAC ###

rule get_bulk_atac:
    output:
        "raw_data/corces_atac/GSE74912_ATACseq_All_Counts.txt",
    shell:
        "sh preprocessing_scripts/download_corces.sh"


rule coembed_10x_atac:
    input:
        "seurat_objects/pbmc_10k_v3.rds",
        "seurat_objects/pbmc_10k_atac.rds"
    output:
        "seurat_objects/pbmc_atac_coembed.rds",
        "analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac/CD4Memory.txt"
    shell:
        "Rscript analysis_code/atac/10x/pbmc_atac_integration.R {input}"

rule pbmc_atac_celltype_coverages:
    input:
        "analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac/CD4Memory.txt"
    output:
        "analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac/CD4Memory.bam"
    conda:
        "envs/atac.yaml"
    shell:
        "sh analysis_code/atac/10x/split_clusters.sh"

rule pbmc_atac_bulk_peak_coverage:
    input:
        "analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac/CD4Memory.bam",
        "raw_data/corces_atac/peaks.bed"
    output:
        "analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac/celltype_coverage.bed"
    shell:
        "sh analysis_code/atac/10x/count_reads_in_peaks.sh"

rule scatac_vs_bulk:
    input:
        "raw_data/corces_atac/GSE74912_ATACseq_All_Counts.txt",
        "analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac/celltype_coverage.bed",
        "analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac/CD4Memory.bam"
    output:
        "figures/atac/pbmc/bulk_correlation_all.pdf"
    shell:
        "Rscript analysis_code/atac/10x/compare_single_cell_bulk.R"

###################### Misc Rules ###########################

rule download_software:
    input:
        "software/download_software.sh"
    output:
        directory("software/seurat-2.3.3/"),
        directory("software/scanorama/"),
        directory("software/seurat")
    shell:
        "sh {input}"
