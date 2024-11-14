package proxy

import (
	"errors"
	"fmt"
	"log"
	"math"
	"math/cmplx"
	"sort"

	"github.com/mjibson/go-dsp/fft"
)

// an enumeration of averaging types (currently mean and median)
type StackType int

const (
	Mean StackType = iota
	Median
)

// Float2D represents a monochrome image using floats.
type Float2D [][]float64

func (c Float2D) Dims() (rows, cols int) {
	return len(c), len(c[0])
}

// creates an empty (zero-valued image)
func NewFloat2D(w, h int) *Float2D {
	f := Float2D(make([][]float64, h))
	for i := 0; i < h; i++ {
		f[i] = make([]float64, w)
	}
	return &f
}

// creates a single-valued image
func NewFloat2DUniform(w, h int, v float64) *Float2D {
	f := Float2D(make([][]float64, h))
	for j := 0; j < h; j++ {
		f[j] = make([]float64, w)
		for i := 0; i < w; i++ {
			f[j][i] = v
		}
	}
	return &f
}

// gets the width and height of the flat array
func (f Float2D) Bounds() (int, int) {
	return len(f), len(f[0])
}

// converts a Float2D array to a list of uint8 (eg, for use as pixel display)
func (f Float2D) AsUint8(rescale bool) ([]uint8, error) {

	H := len(f)
	W := len(f[0])

	if W == 0 || H == 0 {
		return nil, errors.New("array has either zero rows or zero columns")
	}
	allz := make([]float64, W*H)
	px := make([]uint8, W*H)
	mn := 1.0e+10
	mx := -1.0e+10
	for j := 0; j < H; j++ {
		for i := 0; i < W; i++ {
			if f[j][i] > mx {
				mx = f[j][i]
			}
			if f[j][i] < mn {
				mn = f[j][i]
			}
			allz[j*W+i] = f[j][i]
		}
	}
	// we have the range and the values
	// if rescaling, scale the values, then take the uint8
	// otherwise, just take the uint8
	if rescale {
		for i := 0; i < W*H; i++ {
			px[i] = uint8((allz[i] - mn) / (mx - mn) * 255)
		}
	} else {
		for i := 0; i < W*H; i++ {
			px[i] = uint8(allz[i])
		}
	}
	return px, nil
}

// converts a Float2D array to a list of uint16
// these cannot be used for pixel display, as Go image pixels are always stored as uint8
// (except when they aren't)
func (f Float2D) AsUint16(rescale bool) ([]uint16, error) {

	H := len(f)
	W := len(f[0])
	if W == 0 || H == 0 {
		return nil, errors.New("array has either zero rows or zero columns")
	}
	allz := make([]float64, W*H)
	px := make([]uint16, W*H)
	mn := 1.0e+10
	mx := -1.0e+10
	for j := 0; j < H; j++ {
		for i := 0; i < W; i++ {
			if f[j][i] > mx {
				mx = f[j][i]
			}
			if f[j][i] < mn {
				mn = f[j][i]
			}
			allz[j*W+i] = f[j][i]
		}
	}
	// we have the range and the values
	// if rescaling, scale the values, then take the uint8
	// otherwise, just take the uint16

	if rescale {
		for i := 0; i < W*H; i++ {
			px[i] = uint16((allz[i] - mn) / (mx - mn) * 65535)
		}
	} else {
		for i := 0; i < W*H; i++ {
			px[i] = uint16(allz[i])
		}
	}
	return px, nil
}

// shifts the float image by half in x and y, to central DC in amplitude spectra etc.
func (floatIn Float2D) Shift() Float2D {
	// ShowCurrentOperation("Shift")

	w := len(floatIn[0])
	h := len(floatIn)

	var floatOut [][]float64 = make([][]float64, h)
	for j := 0; j < h; j++ {
		floatOut[j] = make([]float64, w)
	}
	for j := 0; j < h; j++ {
		for i := 0; i < w; i++ {
			J := (j + h/2) % h
			I := (i + w/2) % w
			floatOut[J][I] = floatIn[j][i]
		}
	}
	return floatOut
}

// calculates the mean and subtracts from all values
//
// returns detrended data, as well as the mean value
func (f Float2D) RemoveMean() (Float2D, float64) {
	mean := 0.0
	n := float64(len(f) * len(f[0]))
	for j := 0; j < len(f); j++ {
		for i := 0; i < len(f[0]); i++ {
			mean += f[j][i]
		}
	}
	out := make([][]float64, len(f))
	mean /= n
	for j := 0; j < len(f); j++ {
		out[j] = make([]float64, len(f[0]))
		for i := 0; i < len(f[0]); i++ {
			out[j][i] = f[j][i] - mean
		}
	}
	return out, mean
}

// gets the average value
func (f Float2D) Mean() float64 {
	sum := 0.0
	h := len(f)
	w := len(f[0])
	N := float64(w * h)
	for j := 0; j < h; j++ {
		for i := 0; i < w; i++ {
			sum += f[j][i]
		}
	}
	return sum / N
}

// Retains the largest inscribed circle of data centred on the middle of the rectangular array. Other values are
// replaced by either zero or the mean
//
//	usemean - area surrounding circle is replaced by the mean, rather than zero
func (f Float2D) Circularise(usemean bool) Float2D {
	h := len(f)
	w := len(f[0])
	arr := make([][]float64, h)
	// for j := 0; j < h; j++ {
	// 	arr[j] = make([]float64, w)
	// }
	k := h
	if w < k {
		k = w
	}
	var R int
	K := k * k / 4
	// count := 0
	mean := 0.0

	minv := math.MaxFloat64
	for j := 0; j < h; j++ {
		for i := 0; i < w; i++ {
			if f[j][i] < minv {
				minv = f[j][i]
			}
			mean += f[j][j]
		}
	}
	mean /= float64(w * h)
	for j := 0; j < h; j++ {
		arr[j] = make([]float64, w)
		for i := 0; i < w; i++ {
			R = (j-h/2)*(j-h/2) + (i-w/2)*(i-w/2)
			arr[j][i] = f[j][i]
			if R > K {
				if usemean {
					arr[j][i] = mean
				}
				// count++
			}
		}
	}
	return arr
}

// converts linear values to dB values (20*log)
//
//	floor - the minimum value in dB, to prevent large negative values
func (f Float2D) AsDB(floor float64) Float2D {

	h := len(f)
	w := len(f[0])
	if h*w == 0 {
		return nil
	}
	mx := -math.MaxFloat64
	for j := 0; j < h; j++ {
		for i := 0; i < w; i++ {
			if f[j][i] > mx {
				mx = f[j][i]
			}
		}
	}
	out := make([][]float64, h)
	for j := 0; j < h; j++ {
		out[j] = make([]float64, w)
		for i := 0; i < w; i++ {
			out[j][i] = 20 * math.Log10(f[j][i]/mx)
			if out[j][i] > floor {
				out[j][j] = floor
			}
		}
	}
	return out
}

// returns the range of values as minimum, maximum
func (f Float2D) Range() (float64, float64) {
	mx := -math.MaxFloat64
	mn := math.MaxFloat64
	for j := 0; j < len(f); j++ {
		for i := 0; i < len(f[0]); i++ {
			if f[j][i] > mx {
				mx = f[j][i]
			}
			if f[j][i] < mn {
				mn = f[j][i]
			}
		}
	}
	return mn, mx
}

// ranges the data from 0 to 1, then performs log stretch on the ranged numbers, ie, log2(1+k)
func (f Float2D) LogStretch() Float2D {
	h := len(f)
	w := len(f[0])
	var K float64

	mn, mx := f.Range()
	delta := mx - mn
	grey := make([][]float64, h)
	for j := 0; j < h; j++ {
		grey[j] = make([]float64, w)
		for i := 0; i < w; i++ {
			K = (f[j][i] - mn) / (delta)
			// K is between 0 and 1
			K = math.Log2(1 + K)
			// K is between log2(1) and log2(2), ie 0 and 1

			grey[j][i] = K
		}
	}
	return grey
}

// Raises the array to the given power
//
//	k - output A -> A^k
func (f Float2D) Power(k float64) Float2D {
	h := len(f)
	w := len(f[0])

	out := make([][]float64, h)
	for j := 0; j < h; j++ {
		out[j] = make([]float64, w)
		for i := 0; i < w; i++ {
			out[j][i] = math.Pow(f[j][i], k)
		}
	}
	return out
}

// Exponentiates the array to the given power
//
//	k - output A -> e^kA
func (f Float2D) Exponentiate(k float64) Float2D {
	h := len(f)
	w := len(f[0])

	grey := make([][]float64, h)
	for j := 0; j < h; j++ {
		grey[j] = make([]float64, w)
		for i := 0; i < w; i++ {
			grey[j][i] = math.Exp(k * f[j][i])
		}
	}
	return grey
}

// The Normalised Cross-Power Spectrum is useful for detecting the shift between two image arrays.
//
// It results in a phase-only correlation of the images in the frequency domain, which transforms
// to a peak at the maximum lag of the cross-correlation of the two images. This can be picked to
// sub-pixel accuracy, as it has a very strong peak lag which is relatively unaffected by any of
// the neighbouring lags.
func (f Float2D) NormCrossPower(g Float2D) Complex2D {
	F := fft.FFT2Real(f)
	G := fft.FFT2Real(g)
	H := make([][]complex128, len(G))

	for j := 0; j < len(G); j++ {
		H[j] = make([]complex128, len(G[0]))
		for i := 0; i < len(G[0]); i++ {
			H[j][i] = (G[j][i] * cmplx.Conj(F[j][i])) / (F[j][i] * cmplx.Conj(F[j][i]))
		}
	}

	return fft.IFFT2(H)
}

// Linearly interpolates a value at a sub-grid size node position
func (f Float2D) InterpolateBilinear(X, Y float64) float64 {
	x0 := int(X)
	x1 := x0 + 1
	y0 := int(Y)
	y1 := y0 + 1
	dx := X - float64(x0)
	dy := Y - float64(y0)
	a := f[y0][x0]
	b := f[y0][x1]
	c := f[y1][x0]
	d := f[y1][x1]
	E := a + (b-a)*dx
	F := c + (d-c)*dx
	G := E + (F-E)*dy
	return G
}

// Rotates a Float2D by the given angle
//
//	phi - angle by which to rotate the array
func (floatIn Float2D) Rotate(phi float64) *Float2D {
	var x, y, x0, y0, X, Y float64
	var W, H int
	var A float64
	var missedCount int
	var outOfRangeCount int

	w := len(floatIn[0])
	h := len(floatIn)
	c := math.Cos(phi)
	s := math.Sin(phi)
	x0 = float64(w) / 2
	y0 = float64(h) / 2
	diag := math.Sqrt(float64(w*w + h*h))
	W = int(diag)
	H = int(diag)
	X0 := diag / 2
	Y0 := diag / 2

	floatOut := Float2D(make([][]float64, H))
	for j := 0; j < H; j++ {
		floatOut[j] = make([]float64, W)
	}

	for j := 0; j < H; j++ {
		for i := 0; i < W; i++ {
			x = float64(i) - X0
			y = float64(j) - Y0
			X = x*c - y*s + x0
			Y = x*s + y*c + y0
			// only interpolate if the source point is inside the original image
			if X > 0 && X < float64(w)-1 && Y > 0 && Y < float64(h)-1 {

				A = floatIn.InterpolateBilinear(X, Y)

				if A < 0 || A >= 256 {
					outOfRangeCount++
				}
				floatOut[j][i] = A
			} else {
				missedCount++
			}
		}
	}
	fmt.Printf("%d out of range. %d out of rectangle\n", outOfRangeCount, missedCount)
	return &floatOut
}

// shifts an image by an integer number of pixels in x and y
//
//	x,y - shifts to apply
func (f Float2D) TranslatePixel(dx, dy int) *Float2D {
	m := len(f)
	n := len(f[0])
	out := Float2D(make([][]float64, m))
	var I, J int
	for j := 0; j < m; j++ {
		J = j - dy
		out[j] = make([]float64, n)
		for i := 0; i < n; i++ {
			I = i - dx
			if I > -1 && J > -1 && I < n && J < m {
				out[j][i] = f[J][I]
			}
		}
	}
	return &out
}

// shifts an image by a sub-integer number of pixels in x and y using a linear phase shift
//
//	dx,dy - shifts to apply
func (f Float2D) TranslateSubpixelFrequencyDomain(dx, dy float64) Float2D {
	spec := Complex2D(fft.FFT2Real(f))
	return spec.TranslateFD(dx, dy).IFFT().AsReal()
}

func (fp *Float2D) TranslateSubpixelSpaceDomainBilinear(dx, dy float64) Float2D {
	f := *fp
	h := len(f)
	w := len(f[0])
	g := *NewFloat2D(w, h)
	var x, y float64
	for j := 0; j < h; j++ {
		for i := 0; i < w; i++ {
			x = float64(i) - dx
			y = float64(j) - dy
			if x >= 0 && x <= float64(w)-1 && y >= 0 && y <= float64(h)-1 {
				g[j][i] = f.InterpolateBilinear(x, y)
			}
		}
	}
	return g
}

// Aligns two images and report the relative shift by computing the shift from this image to the parameter image.
// Internally, this is done by iteration and can take a while. Use smaller images where posible.
//
//	accuracy - subpixel match precision - iteration stops at this point
//
// returns the shift FROM the base TO the reference
//
//   - The images MUST be equal in size or this will CRASH!
//   - No more than 5 iterations are run, to prevent long run times
func (refp *Float2D) Align(sliderp *Float2D, accuracy float64, maxiterations int) (x, y float64, shifted Float2D, err error) {
	ref := *refp
	slider := *sliderp
	badsize := false
	if len(ref) != len(slider) {
		badsize = true
	}
	if len(ref[0]) != len(slider[0]) {
		badsize = true
	}
	if badsize {
		return 0, 0, nil, errors.New("images are not the same size")
	}
	poc := ref.NormCrossPower(slider)
	dx, dy := poc.MaxPeakFromGaussianFit()
	DX := dx
	DY := dy
	moved := slider.TranslateSubpixelFrequencyDomain(-DX, -DY)
	r := dx*dx + dy*dy
	iterations := 0
	for r > accuracy && iterations < maxiterations {
		iterations++
		poc = ref.NormCrossPower(moved)
		dx, dy := poc.MaxPeakFromGaussianFit()
		r = dx*dx + dy*dy
		DX += dx
		DY += dy
		log.Printf("update %d dx=%03.3f dy=%03.3f cf %03.3f r=%03.3f  DX=%03.3f DY=%03.3f block size %d x %d\n", iterations, dx, dy, accuracy, r, DX, DY, len(slider), len(slider[0]))
		moved = slider.TranslateSubpixelFrequencyDomain(-DX, -DY)
	}
	return DX, DY, moved, nil
}

// returns the power of 2 that is equal to or greater than a number, eg, Get2Power(17)=32
func Get2Power(N int) int {
	n := math.Log2(float64(N))
	padding := math.Abs(n-math.Floor(n)) > 1e-12
	if padding {
		return int(math.Pow(2, math.Floor(n)))
	}
	return N
}

// estimates normalised local variance in a square window, at each pixel of an image
// full windows are carried to the edges for simplicity
//
//	in the returned array, values range from 0,0 to 1.0 (0 means flat, 1 means highly-varying)
func (fp *Float2D) WeightLaplacian() *Float2D {
	f := *fp
	halfwidth := 1
	laplacian := Float2D(make([][]float64, len(f)))
	for j := 0; j < len(f); j++ {
		laplacian[j] = make([]float64, len(f[0]))
	}
	var temp, sum float64
	for j := halfwidth; j < len(f)-halfwidth; j++ {
		for i := halfwidth; i < len(f[0])-halfwidth; i++ {
			sum = 0
			for k := -halfwidth; k < halfwidth; k++ {
				for l := -halfwidth; l < halfwidth; l++ {
					temp = f[j+k][i+l]
					sum += temp
				}
			}
			laplacian[j][i] = sum - 9*f[j][i]
		}
	}
	low, high := laplacian.Range()
	for j := 0; j < len(laplacian); j++ {
		for i := 0; i < len(laplacian[0]); i++ {
			laplacian[j][i] = (laplacian[j][i] - low) / (high - low)
		}
	}
	return &(laplacian)
}

//	*** CHECK THIS CODE ***
//
// calculates local luminance suitability at each pixel of an image
// full windows are carried to the edges for simplicity
//
//	radius - width/height of the square window
func (fp *Float2D) WeightLuminance(radius int) *Float2D {
	f := *fp
	low, high := f.Range()
	log.Printf("Image range for luminance weighting: %1.3f,%1.3f\n", low, high)
	var meanvalue float64
	h := radius / 2
	N := (2*radius + 1)
	N *= N
	lumisuit := Float2D(make([][]float64, len(f)))
	for j := 0; j < len(f); j++ {
		lumisuit[j] = make([]float64, len(f[0]))
	}
	var temp float64
	// σ:=0.2
	// s:=1/σ/σ
	for j := h; j < len(f)-h; j++ {
		for i := h; i < len(f[0])-h; i++ {
			meanvalue = 0
			for k := -h; k < h; k++ {
				for l := -h; l < h; l++ {
					meanvalue += f[j+k][i+l]
				}
			}
			meanvalue /= float64(N * N)
			temp = math.Abs(meanvalue-high/2) / high * 2
			temp *= temp
			temp = math.Exp(-temp * 12.5)
			lumisuit[j][i] = temp
		}
	}
	low, high = lumisuit.Range()

	for j := 0; j < len(lumisuit); j++ {
		for i := 0; i < len(lumisuit[0]); i++ {
			lumisuit[j][i] = (lumisuit[j][i] - low) / high
		}
	}
	// low, high = lumisuit.Range()
	// log.Printf("Range of derived weights: %1.3f,%1.3f\n", low, high)

	return &(lumisuit)
}

// returns the average (mean or median) of supplied arrays, which need to be of the same size.
//
//	method - StackType.Mean or StackType.Median
//	progresschannel - for long-running progresses, percentage completion can be sent to this channel if it exists (connect it to a progress bar, etc)
//
// returns an error if there are no arrays or they are too small (at least 4 x 4 elements)
func AverageFloatArrays(unstacked []Float2D, method StackType, progresschannel *chan float64) (*Float2D, error) {
	N := len(unstacked)
	if N == 0 {
		// WTF? why stack no images?
		return NewFloat2D(1, 1), errors.New("no images are in the slice")
	}

	// if there's only 1 image in the slice, return it
	if N == 1 {
		return &unstacked[0], nil
	}

	// need same dimensions
	var minx int = 1e+10
	var miny int = 1e+10

	base := unstacked[0]
	miny = len(base)
	minx = len(base[0])

	for i := 1; i < len(unstacked); i++ {
		if len(unstacked[i]) < miny {
			miny = len(unstacked[i])
		}
		if len(unstacked[i][0]) < minx {
			minx = len(unstacked[i][0])
		}
	}
	if minx < 4 || miny < 4 {
		return NewFloat2D(1, 1), errors.New("images are too small - need at least 4 pixels in x and y (for arbitrary reasons)")
	}
	var n, sum float64
	n = float64(N)
	// stacking possible
	out := *NewFloat2D(minx, miny)

	for i := 0; i < len(unstacked); i++ {
		if len(unstacked[i]) != miny || len(unstacked[i][0]) != minx {
			return nil, errors.New("arrays must all be the size in both dimensions")
		}
	}

	switch method {
	case Mean:
		for j := 0; j < miny; j++ {
			for i := 0; i < minx; i++ {
				sum = 0
				for k := 0; k < N; k++ {
					sum += unstacked[k][j][i]
				}
				out[j][i] = sum / n
			}
		}
	case Median:
		middle := N / 2
		var holder []float64
		progress := 0.0
		for j := 0; j < miny; j++ {
			progress = float64(j) / float64(miny)
			*progresschannel <- progress
			for i := 0; i < minx; i++ {
				holder = make([]float64, N)
				for k := 0; k < N; k++ {
					// sum += unstacked[k][j][i]
					holder[k] = unstacked[k][j][i]
				}
				out[j][i] = sort.Float64Slice(holder)[middle]
				// out[j][i] = sum / n
			}
		}
	}
	return &out, nil

}

// Converts a 2D Float array to a displayable image.Image Proxy (proxy just allows for metadata to be added).
//
//	bits - 8 or 16, output resolution. Anything else returns an empty image.
func (f Float2D) AsImage(bits int, metadata string, rescale bool) ImageProxy {
	ip := ImageProxy{}
	ip.LoadFromFloats(f, bits, rescale)
	ip.AddMetadata(metadata)
	return ip
}

func (floatIn Float2D) GaussianDownsample() Float2D {

	H := len(floatIn)
	W := len(floatIn[0])
	log.Println(W, H)
	var floatOut, tempRows Float2D
	floatOut = make([][]float64, H/2)
	for j := 0; j < H/2; j++ {
		floatOut[j] = make([]float64, W/2)
	}
	tempRows = make([][]float64, H)
	// tempCols := make([][]float64, W/2)

	for j := 0; j < H; j++ {
		tempRows[j] = make([]float64, W/2)
	}

	// horizontal downsample
	for j := 0; j < H; j++ {
		tempRows[j] = downsample(floatIn[j])
	}

	// vertical downsample
	for i := 0; i < W/2; i++ {
		tempcol := make([]float64, H)
		for j := 0; j < H; j++ {
			tempcol[j] = tempRows[j][i]
		}
		halfcol := downsample(tempcol)

		for j := 0; j < H/2; j++ {
			floatOut[j][i] = halfcol[j]
		}
	}
	return floatOut
}

// extract a rectangle w x h offset from the centre by x,y
func (f Float2D) Middle(w, h, x, y int) Float2D {
	H := len(f)
	W := len(f[0])
	g := make([][]float64, h)
	dx := x - w/2 + W/2
	dy := y - h/2 + H/2

	for j := 0; j < h; j++ {

		g[j] = make([]float64, w)
		for i := 0; i < w; i++ {
			defer func() {
				if r := recover(); r != nil {
					log.Println(w, h, W, H, x, y, i, j, dx, dy)
				}
			}()
			g[j][i] = f[dy+j][dx+i]
		}
	}
	return g
}

var gausscoeffs []float64 = []float64{0.029, 0.235, 0.471, 0.235, 0.029}

func downsample(f []float64) []float64 {
	N := len(f)
	g := make([]float64, N)
	g[0] = f[0]
	g[1] = f[1]
	g[N-1] = f[N-1]
	g[N-2] = f[N-2]
	h := make([]float64, N/2)
	for i := 2; i < N-2; i++ {
		for j := -2; j < 3; j++ {
			g[i] += f[i+j] * gausscoeffs[j+2]
		}
	}

	for i := 0; i < N/2; i++ {
		h[i] = g[2*i]
	}
	return h
}

func (c Float2D) MultiplyElements(a Float2D) (Float2D, error) {
	R, C := c.Dims()
	R2, C2 := a.Dims()
	if R != R2 || C != C2 {
		return nil, errors.New("spectral dimensions do not match")
	}
	var out Float2D = make([][]float64, R)
	for i := 0; i < R; i++ {
		out[i] = make([]float64, C)
		for j := 0; j < C; j++ {
			out[i][j] = c[i][j] * a[i][j]
		}
	}
	return out, nil
}

func (f Float2D) BilateralFilter(sd, sr float64, opradius int) Float2D {
	H := len(f)
	W := len(f[0])

	var out Float2D = make([][]float64, H)

	for j := 0; j < H; j++ {
		out[j] = make([]float64, W)
	}
	log.Println("Made bilateral filter array:",len(out), len(out[0]), W)

	for j := opradius; j < H-opradius-1; j++ {
		for i := opradius; i < W-opradius-1; i++ {
			var divisor float64
			// iterate over the operators
			for k := -opradius; k < opradius+1; k++ {
				for l := -opradius; l < opradius+1; l++ {
					// calculate the weight for this pixel
					gaussian := math.Exp(-(float64(k*k+l*l) / (2 * sd * sd)))
					rangeWeight := math.Exp(-((f[j][i] - f[j+k][i+l]) * (f[j][i] - f[j+k][i+l])) / (2 * sr * sr))
					weight := gaussian * rangeWeight
					out[j][i] += f[j+k][i+l] * weight
					divisor += weight
				}
			}
			out[j][i] /= divisor
		}
	}
	return out
}
