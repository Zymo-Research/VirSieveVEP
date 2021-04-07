# VirSieve VEP

## Summary

This container is part of the Environmental Viral Detection pipeline and covers variant analysis and filtering with GATK.  This operation requires inputs from the _GATK_ portion of this pipeline (specifically the VCF files).

**SECURITY CONCERN**: This pipeline is currently using os.system to run commands and sanitization was causing runs to fail. If running files from untrusted sources, please be sure to sanitize file names to prevent potential command injections into the container.

### File naming and structure
Like other containers in the VirSieve Pipeline, this container is expected to run within a working folder.  This pipeline requires one folder of VCF files for variant annotation with an optional second folder of VCF files that have undergone higher-stringency filtering.  The expected folder name for the standard-filtered VCF files is __filteredVCF__ and for the stringent-filtered VCF files it will be __alignmentArtifactFilteredVCF__.  VCF files will be matched between the folders by their sample name, which is considered the portion of the file name before any dots.
The emphasis of this portion of the pipeline (and the pipeline itself) is to not only identify variants observed in the sample and annotate them, but also to assign them confidence rankings to help filter high-confidence varaints from what is more likely to be biological, chemical, or technical noise in the sample.  Identified variants will be annotated for functional consequences, relative abundances, depths of coverage at their respective loci, and will be assigned a confidence tier between 1 and 3 as follows:

1. High-confidence variants observed with sufficient frequency and passing the high-stringency filter (which is for alignment artifacts by default)
2. Variants of decent confidence that did not pass the high-stringency filter.  Many of these should still be valid.
3. Low-confidence variants that were either poorly supported by the reads (often due to low abundance or low absolute supporting read count) as well as those flagged as having evidence supporting the variant being due to some form of noise (sequencer artifacts, chemical damage pattern, etc.) rather than true biological variation.

### Running the container
To run this container (presumed to be named _virsievevep_ here), simply use the following command:
```bash
docker container run --rm -v /path/to/working/folder:/data virsievevep
```

### Setting non-default options
Some options can be set to non-default values by passing them into the container as environmental variables using the standard Docker commandline technique for setting environmental variables as follows:

| Variable        | Type           | Default  | Description |
| --------------- |:--------------:|:--------:|-------------|
WORKINGFOLDER | string | /data | Working folder name within the container
INPUTFOLDER | string | /$WORKINGFOLDER/filteredVCF | The name of the incoming standard-filtered VCF folder within the working folder
STRINGENTVCFFOLDER | string | /$WORKINGFOLDER/alignmentArtifactFilteredVCF | A folder containing the stringent-filtered VCF files (this file should only have the highest-confidence variants listed)
VEPINTERMEDIATESFOLDER | string | /$WORKINGFOLDER/vepOutputs | Folder for the raw VEP outputs
RESULTSFOLDER | string | /$WORKINGFOLDER/results | Folder for the final outputs to be written


## Contributing

We welcome and encourage contributions to this project from the microbiomics community and will happily accept and acknowledge input (and possibly provide some free kits as a thank you).  We aim to provide a positive and inclusive environment for contributors that is free of any harassment or excessively harsh criticism. Our Golden Rule: *Treat others as you would like to be treated*.

## Versioning

We use a modification of [Semantic Versioning](https://semvar.org) to identify our releases.

Release identifiers will be *major.minor.patch*

Major release: Newly required parameter or other change that is not entirely backwards compatible
Minor release: New optional parameter
Patch release: No changes to parameters

## Authors

- **Michael M. Weinstein** - *Project Lead, Programming and Design* - [michael-weinstein](https://github.com/michael-weinstein)


See also the list of [contributors](https://github.com/Zymo-Research/figaro/contributors) who participated in this project.

## License

This project is licensed under the GNU GPLv3 License - see the [LICENSE](LICENSE) file for details.
This license restricts the usage of this application for non-open sourced systems. Please contact the authors for questions related to relicensing of this software in non-open sourced systems.

## Acknowledgments

We would like to thank the following, without whom this would not have happened:
* The Python Foundation
* The staff at Zymo Research
* The scientific and public health COVID response community
* Our customers


