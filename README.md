![logo](logos/excelerate.png)

# Single cell RNA-seq data analysis with R

This international hands-on course covers several aspects of single cell RNA-seq data analysis, ranging from clustering and differential gene expression analysis to trajectories, cell type identification and spatial transcriptomics. The course is kindly sponsored by the ELIXIR EXCELERATE project.

### [Schedule](schedule.md)

### Date
27.05.2019 - 29.05.2019

### Location
The course is organised in the training room Dogmi at CSC. When you come to the main entrance, turn right to the reception and follow the signs to the course. The street address is Keilaranta 14, Espoo, Finland. You can reach us easily by public transport, please find more details [here](https://www.csc.fi/how-to-reach-us).

### Course computing environment
The software and data required for the exercises have been installed on the classroom computers, and they are also available in CSC's cPouta cloud as a virtual machine image. Please read the [instructions](computing_environment_instructions.md) and how to install the environment on your own computer after the course with [Conda](conda_instructions.md).

## Programme
### Monday 27.5.2019

- [Introduction and experimental design](session-qc/introduction_Jules_GILET.pdf) (Jules Gilet) [Video](https://www.youtube.com/watch?v=BfxDfL1GBzk&list=PLjiXAZO27elC_xnk7gVNM85I2IQl5BEJN&index=2&t=0s)
- [Quality control](session-qc/scRNAseq_QC_Asa_Bjorklund_2019.pdf) (Åsa Björklund) [Video](https://www.youtube.com/watch?v=rOm6UIPhHnc&list=PLjiXAZO27elC_xnk7gVNM85I2IQl5BEJN&index=3&t=0s), [Labs](session-qc/Quality_control.md), [Bonus exercise: Sequencing QC and Ikura pipeline](session-seqmap/sequencing_qc.md)
- [Normalisation](session-normalization/Normalization.pdf) and [removal of confounding factors](session-normalization/confounding-factors.pdf) (Heli Pessa  & Bishwa Ghimire) [Video](https://www.youtube.com/watch?v=gbIks6bA8nI&list=PLjiXAZO27elC_xnk7gVNM85I2IQl5BEJN&index=3), [Video](https://www.youtube.com/watch?v=rhuYhD4GwKw&list=PLjiXAZO27elC_xnk7gVNM85I2IQl5BEJN&index=4), [Labs](session-normalization/Normalization.md), [Labs](session-normalization/confounding-factors.md)
- [Data integration](session-integration/Data_Integration.pdf) (Ahmed Mahfouz) [Video](https://www.youtube.com/watch?v=4KwW90RQz-8&list=PLjiXAZO27elC_xnk7gVNM85I2IQl5BEJN&index=5), [Labs](session-integration/Data_Integration.md) 

### Tuesday 28.5.2019
- [Dimensionality reduction (PCA, tSNE and UMAP)](session-dim-reduction/lecture_dimensionality_reduction.pdf) (Paulo Czarnewski) [Video](https://www.youtube.com/watch?v=qcLJ_JO6bn8&list=PLjiXAZO27elC_xnk7gVNM85I2IQl5BEJN&index=6), [Labs](session-integration/Data_Integration.md), [Rmd file](session-integration/Data_Integration.Rmd)
- Lecture: [Clustering](session-clustering/Clustering.pdf) (Ahmed Mahfouz) [Video](https://www.youtube.com/watch?v=Qa6k7RIwltg&list=PLjiXAZO27elC_xnk7gVNM85I2IQl5BEJN&index=7), [Labs](session-clustering/Clustering.md), [Rmd file](session-clustering/Clustering.Rmd)
- [Differential expression](session-de/session-de.html) (Ståle Nygård) [Video](https://www.youtube.com/watch?v=TibtU3m_bMw&list=PLjiXAZO27elC_xnk7gVNM85I2IQl5BEJN&index=8), [Rmd file 1](session-de/session-de-methods.Rmd), [Rmd file 2](session-de/session-de-methods-evaluation.Rmd)

### Wednesday 29.5.2019
- [Cell type identification](session-celltypeid/celltypeidentification-may2019.pptx) (Philip Lijnzaad) [Video](https://www.youtube.com/watch?v=yjehcKOmqi8&list=PLjiXAZO27elC_xnk7gVNM85I2IQl5BEJN&index=9), [Labs](session-celltypeid/celltypeid.md), [Rmd file](session-celltypeid/celltypeid.Rmd)
- [Trajectories/Pseudo-time I](session-trajectories/trajectory_inference_analysis.pdf) (Paulo, Jules) [Video](https://www.youtube.com/watch?v=XmHDexCtjyw&list=PLjiXAZO27elC_xnk7gVNM85I2IQl5BEJN&index=10), [Labs](session-trajectories/session-trajectories.md), [Labs](session-trajectories/session-trajectories.md#part-ii---diffusion-map)
- [Spatial transcriptomics I](session-spatial/Spatial_transcriptomics_Elixir2019.pdf) (Lars Borm) [Video](https://www.youtube.com/watch?v=XmHDexCtjyw&list=PLjiXAZO27elC_xnk7gVNM85I2IQl5BEJN&index=10), [Video](https://www.youtube.com/watch?v=7m6SGp1fE8w&list=PLjiXAZO27elC_xnk7gVNM85I2IQl5BEJN&index=13&t=0s)

### Prerequisities
In order to participate in this course you should have prior experience in using R.

### Learning objectives
After this course you will be able to:
- use a range of bioinformatics tools to analyze single cell RNA-seq data
- discuss a variety of aspects of single cell RNA-seq data analysis
- understand the advantages and limitations of single cell RNA-seq data analysis

### Lecturers
- Åsa Björklund (NBIS, ELIXIR-SE, Sweden)
- Paulo Czarnewski (NBIS, ELIXIR-SE, Sweden)
- Ahmed Mahfouz (LUMC, Netherlands)
- Ståle Nygård (UIO, Norway)
- Jeongbin Park (Charité-Universitätsmedizin Berlin & de.NBI, Germany)
- Lars Borm (Karolinska Institutet, Sweden)
- Jules Gilet (Institut Curie, France)
- Heli Pessa (University of Helsinki, Finland)
- Bishwa Ghimire (FIMM, Finland)
- Philip Lijnzaad (Princess Maxima Center for Pediatric Oncology, Netherlands)

### Additional information
- Content: <eija.korpelainen@csc.fi>
- Practicalities: <event-support@csc.fi>
