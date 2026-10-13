# Common Errors/Warnings in the NP³ MS Workflow

## Installation and setup

Make sure that all system requirements are installed:

> System requirements: Packages tzdata, zip, wget and which must be installed.

Otherwise, an error with the related missing package may be raised.

----------------------------------

## In the Pre-process Step 2


The following error without interruption:

> error self$finalize()? 

Just ignore this message, this is a bug in the R XCMS 1 from the current version used.

--------------------------------------------

The following warning during the input file (.mzXML) reading followed by the error below. Or just the error:

> Warning: Setting LC_MEASUREMENT failed, using "C"
> Error in pwizModule$open(filename): locale::facet::_S_create_c_locale name not valid

Your system is trying to pass a locale that isn't generated or recognized. 
You can force R to use a universal fallback configuration.

**Solution:**

For quick fix, in the terminal run: 
```{ .text .copy }
export LC_ALL=C.UTF-8
```
or

```{ .text .copy }
export LC_ALL="en_US.UTF-8"
```

For permanent fix, add this environment variable to your default configurations, in the terminal run: 

```{ .text .copy }
echo 'export LC_ALL=C.UTF-8' >> ~/.bashrc
source ~/.bashrc
```

---------------------------------------------

The following error before starting pre-processing:

> Error: package or namespace load failed for ‘utils’:
.onLoad failed in loadNamespace() for 'utils', details:
 call: system(paste(which, shQuote(names[i])), intern = TRUE, ignore.stderr = TRUE)
error: error in running command

> Error: package or namespace load failed for ‘stats’:
.onLoad failed in loadNamespace() for 'utils', details:
 call: system(paste(which, shQuote(names[i])), intern = TRUE, ignore.stderr = TRUE)
error: error in running command

> During startup - Warning messages:
1: package ‘utils’ in options("defaultPackages") was not found  
2: package ‘stats’ in options("defaultPackages") was not found  
Error: package or namespace load failed for ‘dplyr’:
.onLoad failed in loadNamespace() for 'utils', details:
 call: system(paste(which, shQuote(names[i])), intern = TRUE, ignore.stderr = TRUE)
error: error in running command
Execution halted

Your unix system is missing the package `which`.

**Solution:**

Install package `which` in the terminal and try again. Google for installation instructions in your specific 
unix system.
