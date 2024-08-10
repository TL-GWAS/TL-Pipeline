# TarGene-Pipeline

## Documentation Access

Please visit [the docs](https://targene.github.io/targene-pipeline/stable/)

## Running the Workflow

Please refer to the main [NextFlow](https://www.nextflow.io/) documentation for general usage. The main point being that, depending on your cluster specifications, you will need to provide a specific `myprofile` configuration file. If you are part of the University of Edinburgh and simply using Eddie, then the `eddie` profile is already defined. Then simply run:

```bash
nextflow run TARGENE/targene-pipeline -profile myprofile -resume
```

## Live Server

```bash
]activate doc
```

```julia
using LiveServer
servedocs()
```