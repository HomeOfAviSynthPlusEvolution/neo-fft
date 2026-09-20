# Checked plane views and geometry

Specification: MEM-001. The view interface is DS2's `span2d::Plane<T>` at the [reference version](README.md#reference-versions); its required admission and ownership rules are defined below.

## Inputs and outputs

Input is a pointer to live sample objects, sample type T, positive width W and height H, a byte row stride S, and the readable/writable storage extent guaranteed by its owner. Output is `span2d::Plane<const T>` or `span2d::Plane<T>`. The host adapter and validation layer to be implemented must validate the descriptor before constructing the view, using the checks below. Do not introduce a separate plane-view wrapper, owning plane view, or duplicate row/subplane API.

`span2d::Plane<T>` is non-owning. Its constructor takes a byte stride; its `stride()` returns an **element** stride. Width, height and internal element stride use int32. `row_ptr`, `row` and `subplane` do not establish bounds or lifetime safety themselves. A view never outlives its host frame, buffer or workspace lease.

## Validation

Before narrowing, allocation or view construction, require:

1. W,H and S/sizeof(T) are positive and representable as int32.
2. S is divisible by sizeof(T), the pointer is naturally aligned for T, and S >= W*sizeof(T).
3. `(H-1)*S + W*sizeof(T)` is checked and representable in both size_t and ptrdiff_t, and the owner guarantees this accessible extent. Do not add S for a nonexistent final padding row.
4. Every requested subplane satisfies x,y >= 0, w,h > 0, x+w <= W and y+h <= H, with checked arithmetic. Validate the offset before pointer arithmetic.
5. The storage contains live objects of T; raw-byte allocation and reinterpretation alone do not establish object lifetime.

Core views do not support negative stride in this phase. Visible samples need only natural alignment, not SIMD alignment. No kernel may read or write row padding, neighboring planes, inactive batch slots or vector-width tails. A padded owned workspace may expose its padding through a separate explicit extent, never by enlarging a host plane view.

Input, output and scratch allocations are disjoint unless an operator explicitly permits exact in-place operation. Partial overlaps are forbidden. `span2d::RestrictPlane<T>` is an optional internal optimization only after non-aliasing has been established; it must not replace the general adapter type.

## Blocks and spectra

Use a Plane for each 2D slice. A block batch descriptor adds only the active count, capacity, slice/block offsets and FFT shape. It does not reproduce Plane's 2D accessors. A future time dimension is an outer collection of checked slices, not a new stride convention.

Block dimensions are in samples of the **current plane**. FFT3D bw/bh and DFTTest sbsize are not divided by chroma subsampling. Source chroma dimensions come from the host. RGB planes have no YUV centering and no chroma subsampling.

Plan-time arithmetic covers block counts, cover sizes, spectra, batch distances, alignment rounding and all buffer byte totals. Avoid the overflowing ceiling idiom `(a+b-1)/b`; use checked division/remainder with nonnegative a and positive b. Each physical Plane must meet int32 limits even when its owner's total allocation uses size_t.

## Reflection domain

For an axis with source length L and source placed at offset d in a cover of length P, phase 1 admits only `d <= L-1` and `P-d-L <= L-1`. Both pads are nonnegative. This conservatively keeps every reflected read inside the original interval and avoids the references' uninitialized large-padding cases. Reject non-admitted selected-plane geometry at creation; do not infer repeated reflection from an unsafe reference read.

Within the admitted domain, for cover coordinate c let j=c-d. The source coordinate is:

```text
r(j,L) = -j          if j < 0
         j           if 0 <= j < L
         2*L-2-j     if j >= L
```

Apply independently in X and Y. This is reflection **excluding** duplicated endpoints. No clamp, wrap or endpoint repetition is allowed. L=1 works only with zero padding. General repeated reflection for tiny images is a later explicit compatibility decision, not part of phase 1.

## Examples and errors

A uint16 plane W=5,H=3,S=16 has element stride 8 and accessible extent 42 bytes. S=11 is invalid; S=8 cannot hold a visible row. A kernel touching the sixth sample of any row violates the contract even if the owner allocated 48 bytes.

For samples `[10,20,30,40]`, d=2 and P=8, the cover is `[30,20,10,20,30,40,30,20]`. d=4 fails the one-reflection domain. A 1x1 plane with zero pad remains `[10]`.

Known descriptor and allocation-size failures are creation errors; changed or invalid requested-frame descriptors are frame errors. Allocation failures return controlled errors without publishing an incomplete frame. Arithmetic overflow must never be handled by wrapping, truncating or retrying with a smaller image.
