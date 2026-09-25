# FFT direction, half-spectra and normalization

[Contents](../README.md)

The forward transform uses a negative exponent without scaling. The inverse uses a positive exponent and divides by transform volume V:

$$X_k=\sum_{j=0}^{L-1}x_j e^{-2\pi i jk/L},\qquad x_j=\frac1L\sum_{k=0}^{L-1}X_k e^{2\pi i jk/L}.$$

Multidimensional transforms apply this convention on each axis; V is the product of transformed axis lengths. A batch axis contains independent blocks and is not part of V.

## Stored coefficients

Only the final spatial X axis of real input is compressed, to K=floor(Bx/2)+1. The 2D layout is `[By,K]`; DFTTest's 3D layout is `[T,B,K]`, with X contiguous. Y and time remain full axes; even T does not add another compressed axis.

DC has zero frequency on every axis. Even X lengths have a Nyquist column, but arbitrary Y/time coefficients in that column are not necessarily real in multiple dimensions. Conjugate symmetry also reverses the other axes. Filtering uses stored complex coefficients directly; power `re²+im²` is not multiplied by a half-spectrum multiplicity of 2.

## Hand calculation

The educational 1D input `[1,0,0,0]` has full spectrum `[1,1,1,1]` and stored half-spectrum `[1,1,1]`; the normalized inverse reproduces the input. Constant `[2,2,2,2]` has spectrum `[8,0,0,0]`. These are transform examples, not claims about public filter geometry accepting arbitrary 1D sizes.

FFT3D transforms T spatial spectra along time, applying temporal inverse factor 1/T and spatial inverse factor 1/(BxBy) once each. DFTTest uses V=T×B² and an analysis window containing 1/√V. Synthesis uses `U*float(h*V)` in tmode=0, with the synthesis window precomputed at creation, and `(U*V)*h` in tmode=1. These are algebraically equivalent but have different binary32 rounding order. That V compensates the algorithm's scale; it is not redundant normalization.

For DFTTest B=1,T=3 with rectangular windows, amplitudes `[1,2,6]` first receive 1/√3. Retaining the window-shaped mean while clearing residuals yields center output 3, not original center 2. See the full [DFTTest chain](../dfttest.md).

## Paths and rounding

The current implementation uses PocketFFT. `opt=1` controls own operators only; FFT dispatch is separate. Fixed-group batched FFTs, specialized codelets and center-slice inverses implement the same mathematical definitions without changing the meaning of time and batch axes. Different arithmetic paths can round differently; exact constant examples do not prove bitwise equality for all spectra. See [KernelInfo](../kernel-info.md) and [precision boundaries](execution-precision.md).
