!Hard limits on Multiscat parameters included everywhere at compile time

      parameter (NZFIXED_MAX = 550)        !max no of z points in fixed pot.
      parameter (NVFCFIXED_MAX = 262144)      !max no of fourier cmpts (fixed, from file)
      parameter (nmax = 16812)    !NOTE: This is set higher in sams repo, this vm eats too much ram!!      diffraction channels
      parameter (mmax = 550)    !         z grid points
      parameter (nfcx = 262144)    !         max number of potential fourier components per z slice
