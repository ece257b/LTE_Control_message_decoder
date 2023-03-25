# LTE-control-message-decoder
## How to run code
Open up ./codes/rnti_decoder.m in matlab. Change the line 7 to the direct path of "lte_test_R2_2frame.dat", or you can directly run `rnti_decoder <direct_path2lte_test_R2_2frame.dat>` in matlab terminal.

## Software requirement
To run this code, you need to install "LTE Toolbox" in matlab. You can get that from "APPS" section on the top and click "Get More Apps" and search for "LTE Toolbox". It should be included in the UCSD matlab license.  

## Run time of code
It takes about 3 seconds to perform cell search and it will start printing out subframe number, RNTIs and allocation bitmaps. It takes about 10 seconds to finish decoding "lte_test_R2_2frame.dat".

## To play around the tool
I have provided the "./codes/generate_iq.m" files to generate simulated IQ samples from matlab. You can change the RNTI and subframe numbers and feed to "./codes/rnti_decoder.m" again to test with.
