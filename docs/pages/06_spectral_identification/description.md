# Spectral Identification

The NP³ MS Workflow performs the spectral identification against two libraries (both offline):

- The In-Silico predicted MS/MS spectrum of Natural Products 
Database (ISDB) from the Universal Natural Products Database (UNPD) using the tremolo tool (for Unix OS and positive 
ion mode only) - Step 6: command [**tremolo**](tremolo.md);
- And all the GNPS2 community LC experimental data - Step 6.1: command [**gnps_library_search**](gnps_library_search.md).

The command **run** will identify the resulting consensus spectra against these two libraries automatically. After the 
spectral identification a curation of the library annotations is performed for each library result and then 
for the best result from both libraries in the [**Final Identification Curation**](final_identification_curation.md) procedure.

Additionally, the user may also search its result online in the GNPS or GNPS2 website and join the retrieved library 
annotations using the command [**gnps_result**](gnps_result.md). The GNPS curation and the final curation will also be performed here using 
the new joined result.

More details in each command or procedure sections.
