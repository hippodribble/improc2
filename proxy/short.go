package proxy

import (
	"math"

	"github.com/mjibson/go-dsp/fft"
)

type Ushort interface {
	uint8 | uint16
}

type Ushortarray[k Ushort] [][]k

func (g Ushortarray[k]) Range() (int, int) {
	h := len(g)
	w := len(g[0])

	if h == 0 || w == 0 {
		return 0, 0
	}
	mx := -1
	mn := math.MaxInt
	for j := 0; j < h; j++ {
		for i := 0; i < w; i++ {
			if int(g[j][i]) > mx {
				mx = int(g[j][i])
			}
			if int(g[j][i]) < mn {
				mn = int(g[j][i])
			}
		}
	}
	return mn, mx
}

func ListToGrid[k Ushort](g []k, w int) Ushortarray[k] {
	h := len(g) / w

	out := make(Ushortarray[k], h)
	for j := 0; j < h; j++ {
		out[j] = make([]k, w)
		for i := 0; i < w; i++ {
			out[j][i] = g[i+j*w]
		}
	}
	return out
}

func (g Ushortarray[k]) ToList() []k {
	h := len(g)
	w := len(g[0])
	list := make([]k, h*w)
	for j := 0; j < h; j++ {
		for i := 0; i < w; i++ {
			list[i+j*w] = g[j][i]
		}
	}
	return list
}

// convert 8 or 16-bit uint data to float64 for further processsing
func (px Ushortarray[k]) AsFloat() Float2D {
	h := len(px)
	w := len(px[0])
	out := make(Float2D, h)
	for j := 0; j < h; j++ {
		out[j] = make([]float64, w)
		for i := 0; i < w; i++ {
			out[j][i] = float64(px[j][i])
		}
	}
	return out
}

func (sh Ushortarray[k]) Spectrum() Complex2D {
	f,_:=sh.AsFloat().RemoveMean()
	return Complex2D(fft.FFT2Real(f))
}

// converts a 8-bit uint array to 16-bit uint array
func Short8As16(in []uint8)[]uint16{
	out:=make([]uint16,len(in)/2)
	for i:=0;i<len(in)/2;i++{
		out[i]=uint16(in[2*i])*256+uint16(in[2*i+1])
	}
	return out
}


func equaliser[k Ushort](Y []k) *[]k {

	pxout := make([]k, len(Y))
	copy(pxout, Y)
	mapper := make(map[int]int)
	// warm up the map
	for i := 0; i < 256; i++ {
		mapper[i] = 0
	}
	N := len(Y)

	// make histogram
	for i := 0; i < N; i++ {
		I := int((Y)[i])
		mapper[I]++
	}
	cdf := make([]float64, 256)
	cdf[0] = float64(mapper[0])
	// var j uint8
	for i := 1; i < 256; i++ {
		cdf[i] = cdf[i-1] + float64(mapper[i])
	}

	for i := 0; i < 256; i++ {
		mapper[i] = int(cdf[i] * 255 / float64(N))
	}

	for i := 0; i < len(Y); i++ {
		index := int((Y)[i])
		pxout[i] = k(mapper[index])
	}
	return &pxout

}
