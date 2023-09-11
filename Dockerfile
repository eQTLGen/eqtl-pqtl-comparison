FROM nfcore/base
LABEL authors="Urmo Võsa" \
      description="Docker image containing tools for pQTL-eQTL comparisons"

COPY environment.yml /
RUN apt-get update && apt install -y libgmp-dev && apt install -y build-essential
RUN conda env create -f environment.yml && conda clean -a
ENV PATH /opt/conda/envs/pqtlvseqtl/bin:$PATH
ENV TAR="/bin/tar"
RUN ln -s /bin/tar /bin/gtar
COPY temp/hyprcoloc-master.zip /
COPY temp/SumTool-master.zip /
COPY temp/TwoSampleMR-master.zip /
COPY temp/IGUtilityPackage_0.2.98.tar.gz /
RUN R -e "remotes::install_local('hyprcoloc-master.zip')"
RUN R -e "remotes::install_local('SumTool-master.zip')"
RUN R -e "remotes::install_local('TwoSampleMR-master.zip', dependencies = TRUE)"
RUN R -e "install.packages('IGUtilityPackage_0.2.98.tar.gz', repos = NULL, type = 'source', dependencies = TRUE)"