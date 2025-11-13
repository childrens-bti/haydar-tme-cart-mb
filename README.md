# Haydar Lab scRNA-seq data analysis


## To reproduce the code in this repository:
This repository contains a docker image and code used to conduct analyses

1. Clone the repository
```
git clone git@github.com:childrens-bti/haydar-scrna.git
```

2. Pull the docker container:
```
docker pull pgc-images.sbgenomics.com/childrens-bti/haydar-scrna:latest
```
NOTE: if running on a Mac with Apple Silicon chip (M1-M4), please add `--platform linux/amd64`; otherwise add `--platform linux/arm64`

3. Start the docker container, from the `haydar-scrna` folder, run:
```
docker run --platform linux/amd64 --name <CONTAINER_NAME> -d -e PASSWORD=ANYTHING -p 8787:8787 -v $PWD:/home/rstudio/haydar-scrna pgc-images.sbgenomics.com/childrens-bti/haydar-scrna:latest
```
NOTE: if running on a Mac with Apple Silicon chip (M1-M4), please add `platform linux/amd64`

4. To execute shell within the docker image, from the `haydar-scrna` folder, run:
```
docker exec -ti <CONTAINER_NAME> bash
```

5. Navigate to an analysis module and run the shell script:
```
cd /home/rstudio/haydar-scrna/analyses/module_of_interest
```

### Below is the level one directory structure listing the analyses and data files used in this repository

```
.
├── analyses
├── data
├── Dockerfile
├── docs
├── figures
├── LICENSE
├── README.md
├── download_data.sh
└── scripts
```

## Code Authors

Sam Chen ([@sychen9584](https://github.com/sychen9584))
