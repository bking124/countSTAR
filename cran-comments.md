## Resubmission
This is a resubmission. In this version I have changed several things according to the comments given by the CRAN team:

* Added several references to the DESCRIPTION
* Added \value to certain .Rd files which were missing this field
* Exported two functions which were previously unexported
* Removed examples for all other unexported functions
* Removed \dontrun for two examples with small compute time
* Changed \dontrun to \donttest for examples that take >5 sec

## R CMD check results

0 errors | 0 warnings | 1 note

* The note seems to regard a URL (to a paper) which I have checked to be valid on all browsers and systems to which I have access.
