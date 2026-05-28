Based on CMSSW\_16\_1\_0\_pre4

## Running the vertexing 

cd PrimaryVertexProducer/test
cmsRun PrimaryVertexProducer/test/testVertexing.py -m [MODE]  -i [INPUT FILES] -o [OUTPUT FILE] 

where
[MODE] can be "Old" "InBlocks" "Alpaka"
[INPUT FILES] can either be a single root file or a .txt with a set

DQM output can be added with the '-d' option

For example
cmsRun testVertexing.py -m Alpaka  -i Zmumu\_16X.txt -o Zmumu\_RelVal\_16X.root -d

## Harvesting DQM plots:

cmsRun harvester.py [INPUT FILE]

## Making DQM plots from harvested files:

makeTrackValidationPlots.py --extended [FILENAME1] [FILENAME2] ... 
