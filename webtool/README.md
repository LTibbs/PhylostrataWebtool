# MaizeGDB Phylostrata Webtool
This subfolder contains the scripts used for the webtool to visualize phylostrata data at MaizeGDB. Use the two json files from the [data_processing](https://github.com/LTibbs/PhylostrataWebtool/tree/main/data_processing) subfolder as input. 

In this code, I have used `NOTE` in comments to denote places where the user may want or need to make changes to customize the results. For example, the user can specify the example species for each phylostratum and provide example images.

In the main folder contains scripts used to make the webtool main pages, including the search box, downloads, etc.

There are also three subfolders:
1. [gene_pages](https://github.com/LTibbs/PhylostrataWebtool/tree/main/webtool/gene_pages) contains the scripts to create the gene details pages and associated BLAST results pages.
2. [image_code](https://github.com/LTibbs/PhylostrataWebtool/tree/main/webtool/image_code) contains the scripts to create the visualization of the phylostrata results.
3. [other_images](https://github.com/LTibbs/PhylostrataWebtool/tree/main/webtool/other_images) contains all of the other files referenced in the code, for example images of the strata example species, results files for the Downloads page, etc. 

## Deployment
The Phylostrata tool has been containerized into a docker image and is available at <a href="https://hub.docker.com/r/maizegdb/phylostrata">docker.io/maizegdb/phylostrata</a>. MaizeGDB is running the tool as a container on our Podman server, but it can also be deployed locally. If you have docker installed on your device, you can download and run the tool with this command in the terminal:

`$ docker run -it -d -p 8080:80 --name phylostrata maizegdb/phylostrata`

Then open a web browser and navigate to http://localhost:8080. You should have full access to the tool from your local device here. To stop the container run this command:

`$ docker stop phylostrata`

The container will stop running with this command, but the image is still downloaded to your device and is taking up disk space. To remove the image, run this command:

`$ docker rm phylostrata`
