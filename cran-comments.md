## gclm 0.0.1

Re-submission to CRAN 

* Added reference to arxiv preprint in DESCRIPTION.
* Added `"cph"` in Authors@R 

There are notes on some checks below:

* Possibly mis-spelled words in DESCRIPTION:
  Lyapunov (3:29, 13:30)
  proximal (15:21)
  
  **Both words are correctly spelled**
  
* checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    'examples_i386' 'examples_x64' 'gclm-Ex_i386.Rout' 'gclm-Ex_x64.Rout'
    'tests_i386' 'tests_x64'

## Test environment

* Ubuntu 18.04.2      (64-bit)  R 4.0.0 (local) 
* Windows Server 2008 (64-bit)  R 4.0.0 (win-builder, r-release)
* Windows Server 2008 (64-bit)  R 3.6.3 (win-builder, r-oldrelease)
* Windows Server 2008 (64-bit)  R devel (win-builder, r-devel)
* Ubuntu Linux 16.04 LTS, R-release, GCC         (R-hub)
* Debian Linux, R-devel, GCC ASAN/UBSAN          (R-hub)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (R-hub)
* Fedora Linux, R-devel, clang, gfortran         (R-hub)

## R CMD Check results

### Ubuntu 18.04.2 R 3.6.3 (local) 

#### devtools::check()

Duration: 14.4s

0 errors | 0 warnings | 0 note

#### R CMD check --as-cran

Status: 1 NOTE

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Gherardo Varando <gherardo.varando@gmail.com>’

New submission

### win-builder

#### R-release

Installation time in seconds: 10
Check time in seconds: 54
Status: 1 NOTE

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Gherardo Varando <gherardo.varando@gmail.com>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  Lyapunov (3:29, 13:30)
  proximal (15:21)

#### R-devel 

Installation time in seconds: 10
Check time in seconds: 54
Status: 1 NOTE

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Gherardo Varando <gherardo.varando@gmail.com>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  Lyapunov (3:29, 13:30)
  proximal (15:21)

#### R-oldrelease

Installation time in seconds: 11
Check time in seconds: 51
Status: 1 NOTE

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Gherardo Varando <gherardo.varando@gmail.com>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  Lyapunov (3:29, 13:30)
  proximal (15:21)

### r-hub

#### Ubuntu Linux 16.04 LTS, R-release, GCC

0 errors  | 0 warnings  | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Gherardo Varando <gherardo.varando@gmail.com>’
  
  New submission
  
  Possibly mis-spelled words in DESCRIPTION:
    Lyapunov (3:29, 13:30)
    proximal (15:21)
    
#### Debian Linux, R-devel, GCC ASAN/UBSAN

0 errors  | 0 warnings  | 0 notes 

#### Windows Server 2008 R2 SP1, R-devel, 32/64 bit

0 errors  | 0 warnings  | 2 notes

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Gherardo Varando <gherardo.varando@gmail.com>'
  
  New submission
  
  Possibly mis-spelled words in DESCRIPTION:
    Lyapunov (3:29, 13:30)
    
* checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    'examples_i386' 'examples_x64' 'gclm-Ex_i386.Rout' 'gclm-Ex_x64.Rout'
    'tests_i386' 'tests_x64'
    
#### Fedora Linux, R-devel, clang, gfortran

0 errors  | 0 warnings  | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Gherardo Varando <gherardo.varando@gmail.com>’
  
  New submission
  
  Possibly mis-spelled words in DESCRIPTION:
    Lyapunov (3:29, 13:30)
