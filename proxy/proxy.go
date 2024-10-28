package proxy

import (
	"errors"
	"fmt"
	"image"
	"image/color"
	"image/draw"
	_ "image/gif"
	_ "image/jpeg"
	_ "image/png"
	"log"
	"math"
	"os"
	"strings"
	"sync"
	"time"

	"github.com/disintegration/imaging"
	"github.com/mjibson/go-dsp/fft"
)

type Component int

const (
	R Component = iota
	G
	B
	A
	Y
	Cb
	Cr
)

// ImageProxy embeds image.Image, adding functions for image processing. Metadata can be appended at each stage to create a processing history
//
//	metadata  add a line of text when the image is modified to track its history
//	Selected  useful for display windows, to check whether this specific ImageProxy will be used
//	Path      the original file path if the ImageProxy was loaded from a resource
type ImageProxy struct {
	image.Image
	Config   image.Config
	metadata []string
	Selected bool
	Path     string
}

// this interface is used to ensure via assertion that an image type supports the Set() method, which is otherwise not part of the image>Image interface
type Changeable interface {
	Set(x, y int, c color.Color)
}

// adds another line of metadata to the ImageProxy
func (ip *ImageProxy) AddMetadata(data string) ImageProxy {
	if ip.metadata == nil {
		ip.metadata = make([]string, 0)
	}
	t := time.Now().Format("2006-02-01 15:04:05")
	ip.metadata = append(ip.metadata, t+"  "+data)
	return *ip
}

// returns all metadata created for the ImageProxy as a single string.
func (ip ImageProxy) AllMetadata() string {
	s := ""
	if ip.metadata == nil {
		return s
	}
	for i := 0; i < len(ip.metadata); i++ {
		s += ip.metadata[i] + "\n"
	}
	return s
}

// Returns the last element of the metadata slice, or an empty string if the slice is empty.
func (ip ImageProxy) LastMetadata() string {
	if ip.metadata == nil {
		return ""
	}
	return ip.metadata[len(ip.metadata)-1]
}

// returns the type of the underlying image.Image
func (ip ImageProxy) Type() string {
	Type := "Unknown"
	switch ip.Image.(type) {
	case *image.Gray:
		Type = "Gray 8 bit"
	case *image.Gray16:
		Type = "Gray 16 bit"
	case *image.RGBA:
		Type = "RGBA 8 bit"
	case *image.NRGBA:
		Type = "NRGBA 8 bit"
	case *image.RGBA64:
		Type = "RGBA 16 bit"
	case *image.NRGBA64:
		Type = "NRGBA 16 bit"
	case *image.YCbCr:
		Type = "YCbCr"
	case *image.NYCbCrA:
		Type = "NYCbCrA"
	case *image.CMYK:
		Type = "CMYK"
	case *image.Paletted:
		Type = "Paletted"
	case *image.Alpha:
		Type = "Alpha"
	case *image.Alpha16:
		Type = "Alpha16"
	}
	return Type
}

// adds a colour palette to an 8-bit grayscale image. If the input is not 8-bit gray, it just returns the input
//
//	p - colour palette (just a []color.Color 256 in length)
func (ip *ImageProxy) GrayToPaletted(p color.Palette) *ImageProxy {
	paletted := image.NewPaletted(ip.Bounds(), p)

	switch im := ip.Image.(type) {
	case *image.Gray:
		paletted.Pix = im.Pix
		op := ImageProxy{Image: paletted}
		op.AddMetadata("Converted from Grayscale")
		return &op
	default:
		return ip
	}
}

// removes the colour palette from an image to make it 8-bit Gray
// - if the image is not an 8-bit paletted image, it returns the otiginal image
func (ip *ImageProxy) PalettedToGray() *ImageProxy {
	switch im := ip.Image.(type) {
	case *image.Paletted:
		newgray := image.NewGray(ip.Bounds())
		newgray.Pix = im.Pix
		op := new(ImageProxy)
		op.Image = newgray
		return op
	default:
		return ip
	}
}

// loads an image from a file
func (ip *ImageProxy) LoadFromFile(path string) error {
	r, err := os.Open(path)
	if err != nil {
		return errors.New("could not open the request URI")
	}
	defer r.Close()

	config, _, err := image.DecodeConfig(r)
	if err != nil {
		return errors.New("decode config: " + err.Error())
	}
	ip.Config = config

	r.Seek(0, 0)

	im, _, err := image.Decode(r)
	if err != nil {
		return errors.New("decode image: " + err.Error())
	}

	ip.Image = im

	// if ip.Type() == "YCbCr" {
	// 	ip.Image=ip.YCbCrToRGBA()
	// 	ip.AddMetadata("Cnverted YCbCr to RGBA")
	// 	log.Println("Transformed Ycbcr to ", ip.Type())
	// 	return nil
	// }
	return nil
}

// generates a grayscale image from the provided floats, autoscaling amplitudes to fit the range available.
// Because the input is a simple array of Float2D, only Gray or Gray16 images are returned.
//
//	f - Float2D array of data.
//	nbits    - desired resolution of image in bits, either 8 or 16. No other options are allowed.
//	rescale  - whether to rescale the data to 0..255 or 0..65535 range - if false, data will be truncated to these ranges for safety
func (ip *ImageProxy) LoadFromFloats(floatIn Float2D, nbits int, rescale bool) error {
	// f := *floats
	if nbits != 8 && nbits != 16 {
		return errors.New("specify either 8 or 16-bit resolution")
	}

	var v float64
	var lo, hi uint8
	upper := uint16(255 << 8)

	h := len(floatIn)
	w := len(floatIn[0])

	if !rescale {
		if nbits == 8 {
			im := image.NewGray(image.Rect(0, 0, w, h))
			for j := 0; j < h; j++ {
				for i := 0; i < w; i++ {
					v = floatIn[j][i]
					if v < 0 {
						v = 0
					}
					if v > 255 {
						v = 255
					}
					im.Pix[j*w+i] = uint8(v)
				}
			}
			ip.Image = im
			ip.Config = image.Config{ColorModel: im.ColorModel(), Height: h, Width: w}
			ip.AddMetadata("Gray 16 created from floating array with no rescaling")
		} else {
			im := image.NewGray16(image.Rect(0, 0, w, h))
			for j := 0; j < h; j++ {
				for i := 0; i < w; i++ {
					v = floatIn[j][i]
					if v < 0 {
						v = 0
					}
					if v > 65535 {
						v = 65535
					}
					hi = uint8((uint16(v) & upper) >> 8)
					lo = uint8(uint16(v) & 255)
					im.Pix[2*(j*w+i)] = hi
					im.Pix[2*(j*w+i)+1] = lo
				}
			}
			ip.Image = im
			ip.Config = image.Config{ColorModel: im.ColorModel(), Height: h, Width: w}
			ip.AddMetadata("Gray 16 created from floating array with no rescaling")
		}
	} else {
		// rescaling all values to the appropriate range (8 or 16 bits)
		mn, mx := floatIn.Range()
		rng := mx - mn
		if nbits == 8 {
			im := image.NewGray(image.Rect(0, 0, w, h))
			for j := 0; j < h; j++ {
				for i := 0; i < w; i++ {
					v = (floatIn[j][i] - mn) / rng * 255
					if v < 0 {
						v = 0
					}
					if v > 255 {
						v = 255
					}
					im.Pix[j*w+i] = uint8(v)
				}
			}
			ip.Image = im
			ip.Config = image.Config{ColorModel: im.ColorModel(), Height: h, Width: w}
			ip.AddMetadata("Gray 16 created from floating array with range rescaling")
		} else {
			im := image.NewGray16(image.Rect(0, 0, w, h))
			for j := 0; j < h; j++ {
				for i := 0; i < w; i++ {
					v = (floatIn[j][i] - mn) / rng * 65535
					if v < 0 {
						v = 0
					}
					if v > 65535 {
						v = 65535
					}
					hi = uint8((uint16(v) & upper) >> 8)
					lo = uint8(uint16(v) & 255)
					im.Pix[2*(j*w+i)] = hi
					im.Pix[2*(j*w+i)+1] = lo
				}
			}
			ip.Image = im
			ip.Config = image.Config{ColorModel: im.ColorModel(), Height: h, Width: w}
			ip.AddMetadata("Gray 16 created from floating array with range rescaling")
		}
	}
	return nil
}

// makes a 32 bit image from r,g,b components (alpha is 255)
//
// rescaling means that ALL components are rescaled with the same factor, to prevent coplour changes
// no rescaling just tries to convert to 8-bit. Everything outside of 0-255 is truncated.
func (ip *ImageProxy) CreateFromRGB(RGB []Float2D, rescale bool) error {

	if len(RGB) < 3 {
		return errors.New("need 3 planes in order to create an image")
	}

	pmax := 0.0
	pmin := 255.0

	w := len(RGB[0][0])
	h := len(RGB[0])
	px := make([]uint8, 4*w*h)
	imout := image.NewNRGBA(image.Rect(0, 0, w, h))

	if !rescale {
		for j := 0; j < len(RGB[0]); j++ {
			for i := 0; i < w; i++ {
				for k := 0; k < len(RGB); k++ {
					if RGB[k][j][i] >= 0 && RGB[k][j][i] <= 255 {
						px[4*(j*w+i)+k] = uint8(RGB[k][j][i])
					}
				}
				px[4*(j*w+i)+3] = 255
			}
		}
		imout.Pix = px
		ip.Image = imout
		return nil
	}

	for i := range RGB {
		plane := RGB[i]
		for j := range plane {
			row := plane[j]
			for k := range row {
				if row[k] > pmax {
					pmax = row[k]
				}
				if row[k] < pmin {
					pmin = row[k]
				}
			}
		}
	}

	var a, b float64
	scope := pmax - pmin

	for j := 0; j < len(RGB[0]); j++ {
		for i := 0; i < w; i++ {
			for k := 0; k < len(RGB); k++ {
				a = RGB[k][j][i]
				b = (a - pmin) / scope * 255
				if b < 0 {
					b = 0
				}
				if b > 255 {
					b = 255
				}
				px[4*(j*w+i)+k] = uint8(b)
			}
			px[4*(j*w+i)+3] = 255
		}
	}
	imout.Pix = px
	ip.Image = imout
	return nil

}

// generates an NRGBA or NRGBA64 image from the provided Float2D arrays. Data can be either clipped or rescaled
func (ip *ImageProxy) LoadFromFloatsAsRGBAComponents(RGBA []Float2D, nbits int, rescale bool) error {
	if len(RGBA) != 4 {
		return errors.New("need precisely 4 Float2D arrays of the same size")
	}
	h := len(RGBA[0])
	w := len(RGBA[0][0])
	for i := 1; i < 4; i++ {
		if len(RGBA[i]) != h || len(RGBA[i][0]) != w {
			return errors.New("arrays are not the same size")
		}
	}
	if nbits == 8 {
		u8c := make([][]uint8, 4)
		for i := 0; i < 4; i++ {
			x, err := RGBA[i].AsUint8(rescale)
			if err != nil {
				return err
			}
			u8c[i] = x
		}
		imout := image.NewNRGBA(image.Rect(0, 0, w, h))
		for j := 0; j < h; j++ {
			for i := 0; i < w; i++ {
				imout.SetNRGBA(i, j, color.NRGBA{u8c[0][j*w+i], u8c[1][j*w+i], u8c[2][j*w+i], u8c[3][j*w+i]})
			}
		}
		ip.Image = imout
		ip.AddMetadata("8-bit NRGBA from float components")

	} else {
		u16c := make([][]uint16, 4)
		for i := 0; i < 4; i++ {
			x, err := RGBA[i].AsUint16(rescale)
			if err != nil {
				return err
			}
			u16c[i] = x
		}
		imout := image.NewNRGBA64(image.Rect(0, 0, w, h))
		for j := 0; j < h; j++ {
			for i := 0; i < w; i++ {
				imout.SetNRGBA64(i, j, color.NRGBA64{u16c[0][j*w+i], u16c[1][j*w+i], u16c[2][j*w+i], u16c[3][j*w+i]})
			}
		}
		ip.Image = imout
		ip.AddMetadata("16-bit RGBA from float components")
	}
	return nil
}

// returns the requested RGBA component as a Float2D array
//
//	c - component, either R,G,B or A (of type Component, an iota)
func (ip ImageProxy) GetComponentAsFloat2(c Component) *Float2D {
	w := ip.Bounds().Dx()
	h := ip.Bounds().Dy()
	var r, g, b, a uint32
	log.Println(w, h, c)

	floatOut := Float2D(make([][]float64, h))

	for j := 0; j < h; j++ {
		floatOut[j] = make([]float64, w)
	}
	switch c {
	case R:
		for j := 0; j < h; j++ {
			for i := 0; i < w; i++ {
				r, _, _, _ = ip.At(i, j).RGBA()
				floatOut[j][i] = float64(r)
			}
		}
	case G:
		for j := 0; j < h; j++ {
			for i := 0; i < w; i++ {
				_, g, _, _ = ip.At(i, j).RGBA()
				floatOut[j][i] = float64(g)
			}
		}
	case B:
		for j := 0; j < h; j++ {
			for i := 0; i < w; i++ {
				_, _, b, _ = ip.At(i, j).RGBA()
				floatOut[j][i] = float64(b)
			}
		}
	case A:
		for j := 0; j < h; j++ {
			for i := 0; i < w; i++ {
				_, _, _, a = ip.At(i, j).RGBA()
				floatOut[j][i] = float64(a)
			}
		}
	}
	return &floatOut
}

// returns the requested RGBA component as a Float2D array
//
//	c - component, either R,G,B or A (of type Component, an iota)
func (ip ImageProxy) GetComponentAsFloat(c Component) *Float2D {
	w := ip.Bounds().Dx()
	h := ip.Bounds().Dy()
	var r, g, b, a uint32
	log.Println(w, h, c)

	floatOut := Float2D(make([][]float64, h))

	for j := 0; j < h; j++ {
		floatOut[j] = make([]float64, w)
	}

	for j := 0; j < h; j++ {
		for i := 0; i < w; i++ {
			r, g, b, a = ip.At(i, j).RGBA()
			switch c {
			case R:
				floatOut[j][i] = float64(r)
			case G:
				floatOut[j][i] = float64(g)
			case B:
				floatOut[j][i] = float64(b)
			case A:
				floatOut[j][i] = float64(a)
			}
		}
	}

	return &floatOut
}

// returns a single component of an image as a grayscale image, either 8 or 16 bits (automatic based on input image)
func (ip *ImageProxy) GetComponent(c Component) *ImageProxy {
	w := ip.Bounds().Dx()
	h := ip.Bounds().Dy()
	var r, g, b, a uint32
	op := new(ImageProxy)
	if strings.HasSuffix(ip.Type(), "16") || strings.HasSuffix(ip.Type(), "64") {
		// 16 bits required
		im := image.NewGray16(image.Rect(0, 0, w, h))
		switch c {
		case R:
			for j := 0; j < h; j++ {
				for i := 0; i < w; i++ {
					r, _, _, _ = ip.At(i, j).RGBA()
					im.Set(i, j, color.Gray16{uint16(r)})
				}
			}
			op.AddMetadata("Extracted red component")
		case G:
			for j := 0; j < h; j++ {
				for i := 0; i < w; i++ {
					_, g, _, _ = ip.At(i, j).RGBA()
					im.Set(i, j, color.Gray16{uint16(g)})
				}
			}
			op.AddMetadata("Extracted green component")
		case B:
			for j := 0; j < h; j++ {
				for i := 0; i < w; i++ {
					_, _, b, _ = ip.At(i, j).RGBA()
					im.Set(i, j, color.Gray16{uint16(b)})
				}
			}
			op.AddMetadata("Extracted blue component")
		case A:
			for j := 0; j < h; j++ {
				for i := 0; i < w; i++ {
					_, _, _, a = ip.At(i, j).RGBA()
					im.Set(i, j, color.Gray16{uint16(a)})
				}
			}
			op.AddMetadata("Extracted alpha component")
		}
		op.Image = im

	} else {
		// 8 bits only
		im := image.NewGray(image.Rect(0, 0, w, h))

		switch c {
		case R:
			for j := 0; j < h; j++ {
				for i := 0; i < w; i++ {
					r, _, _, _ = ip.At(i, j).RGBA()
					im.Pix[j*w+i] = uint8(r)
				}
			}
			op.AddMetadata("Extracted red component")
		case G:
			for j := 0; j < h; j++ {
				for i := 0; i < w; i++ {
					_, g, _, _ = ip.At(i, j).RGBA()
					im.Pix[j*w+i] = uint8(g)
				}
			}
			op.AddMetadata("Extracted green component")
		case B:
			for j := 0; j < h; j++ {
				for i := 0; i < w; i++ {
					_, _, b, _ = ip.At(i, j).RGBA()
					im.Pix[j*w+i] = uint8(b)
				}
			}
			op.AddMetadata("Extracted blue component")
		case A:
			for j := 0; j < h; j++ {
				for i := 0; i < w; i++ {
					_, _, _, a = ip.At(i, j).RGBA()
					im.Pix[j*w+i] = uint8(a)
				}
			}
			op.AddMetadata("Extracted alpha component")
		}
		op.Image = im
	}
	return op
}

// converts to RGBA
func (ip *ImageProxy) YCbCrToRGBA() *image.RGBA {
	// op := new(ImageProxy)
	im := ip.Image
	if iy, ok := im.(*image.YCbCr); ok {
		outim := image.NewRGBA(im.Bounds())
		// 24 bit image

		// Y:=iy.Y
		// Cb:=iy.Cb
		// Cr:=iy.Cr
		// log.Printf("No of elements Y: %d Cb: %d Cr: %d\n",len(Y),len(Cb),len(Cr))
		// var r, g, b, a uint32
		for i := 0; i < ip.Bounds().Dx(); i++ {
			for j := 0; j < ip.Bounds().Dy(); j++ {
				c := iy.YCbCrAt(i, j)
				outim.Set(i, j, c)
			}
		}
		// op.Image = outim
		// log.Println("Successful transform of colours!")
		// log.Println(op.Type(),outim.ColorModel())
		return outim
	}

	return nil
}

// copies an existing image, regardless of underlying type
func (src *ImageProxy) Copy() *ImageProxy {
	bounds := src.Bounds()
	var dst draw.Image
	op := new(ImageProxy)

	switch src.Image.(type) {
	case *image.Gray:
		dst = image.NewGray(bounds)
	case *image.Gray16:
		dst = image.NewGray16(bounds)
	case *image.RGBA:
		dst = image.NewRGBA(bounds)
	case *image.RGBA64:
		dst = image.NewRGBA64(bounds)
	case *image.NRGBA:
		dst = image.NewNRGBA(bounds)
	case *image.NRGBA64:
		dst = image.NewNRGBA64(bounds)
	// case *image.YCbCr:
	// 	dst = image.NewYCbCr(bounds, src.(*image.YCbCr).SubsampleRatio)
	case *image.CMYK:
		dst = image.NewCMYK((bounds))
	default:
		dst = image.NewRGBA(bounds)
	}

	draw.Draw(dst, bounds, src, bounds.Min, draw.Src)
	op.Image = dst
	op.Config = image.Config{ColorModel: dst.ColorModel(), Width: dst.Bounds().Dx(), Height: dst.Bounds().Dy()}
	op.AddMetadata("Copied as " + op.Type())
	return op

}

// masks the outside of a maximal inscribed circle in the image
func (ip *ImageProxy) Circularise() *ImageProxy {
	ipout := ip.Copy()
	w := ip.Bounds().Dx()
	h := ip.Bounds().Dy()
	s := h // size of a side of the square output
	if w < h {
		s = w
	}
	r := s * s / 4

	if zoid, ok := ipout.Image.(Changeable); ok {

		var x2, y2 int
		for j := 0; j < h; j++ {
			y2 = j - h/2
			y2 *= y2
			for i := 0; i < w; i++ {
				x2 = i - w/2
				x2 *= x2
				if x2+y2 > r {
					zoid.Set(i, j, color.Black)
				}
			}
		}

	}
	ipout.Config = image.Config{ColorModel: ipout.Image.ColorModel(), Width: ipout.Bounds().Dx(), Height: ipout.Bounds().Dy()}
	ipout.AddMetadata("Circularised")
	return ipout
}

// converts luminance values in any colour space to a Float2D array
func (ip *ImageProxy) LuminanceAsFloat() *Float2D {
	var r, g, b uint32
	// fd:=ip.Image
	// convert pixels directly to luminance values as uint16

	// fmt.Println(ip.AllMetadata())

	m := Float2D(make([][]float64, ip.Bounds().Dy()))
	for j := 0; j < ip.Bounds().Dy(); j++ {
		m[j] = make([]float64, ip.Bounds().Dx())
		for i := 0; i < ip.Bounds().Dx(); i++ {
			r, g, b, _ = ip.At(i, j).RGBA()
			// if r+g+b!=0{
			// 	fmt.Println(i,j,r,g,b)
			// }
			m[j][i] = float64(max(r, g, b))

		}
	}
	return &m
}

// returns the spectrum of the image's luminance values
func (ip *ImageProxy) FFT() *Complex2D {
	a := Complex2D(fft.FFT2Real(*ip.LuminanceAsFloat()))
	return &a
}

// apply histogram equalisation to the luminance channel of an image
func (ip *ImageProxy) HistogramEqualise() *ImageProxy {

	var mmod ImageProxy

	switch ip.Image.(type) {
	case *image.Gray:
		im := ip.Image.(*image.Gray)
		newim := image.NewGray(ip.Bounds())
		// px := equalise8bit(&im.Pix)
		px := equaliser(im.Pix)
		newim.Pix = *px
		mmod = ImageProxy{Image: newim}
		mmod.AddMetadata("Equalised")

	case *image.Gray16:
		fmt.Println("16-bit")
		im := ip.Image.(*image.Gray16)
		gray16pix := convert8to16(&im.Pix)
		histeq16 := equalise16(gray16pix)
		// px := equalise16bit(&im.Pix)
		newim := image.NewGray(ip.Bounds())
		newim.Pix = *convert16to8(histeq16)
		mmod = ImageProxy{Image: newim}
		mmod.AddMetadata("Equalised - 16-bit")

	case *image.RGBA:
		// cheaty cheaty - take as v the max of r/g/b - store as a histogram, make the cdf and map the values from the array to these.
		im := ip.Image.(*image.RGBA)
		w := im.Bounds().Dx()
		h := im.Bounds().Dy()
		N := w * h
		var r, g, b, a, v uint8
		mapper := make(map[int]int)
		for i := 0; i < N; i++ {
			r = im.Pix[4*i+0]
			g = im.Pix[4*i+1]
			b = im.Pix[4*i+2]
			v = max(r, g, b)
			mapper[int(v)]++
		}
		// our mapper is complete
		cdf := make([]float64, 256)
		cdf[0] = float64(mapper[0])
		for i := 1; i < 256; i++ {
			cdf[i] = float64(mapper[i]) + cdf[i-1]
		}
		// cdf complete
		// replace the mapper with the cdf function values
		for i := 0; i < 256; i++ {
			mapper[i] = int(cdf[i] * 255 / float64(N))
		}
		imout := image.NewRGBA(im.Bounds())

		for i := 0; i < N; i++ {

			r = im.Pix[4*i+0]
			g = im.Pix[4*i+1]
			b = im.Pix[4*i+2]
			a = im.Pix[4*i+3]
			v1 := float64(max(r, g, b))
			// v is mapped to v2 by the equalisation
			v2 := float64(mapper[int(v1)])
			k := v2 / v1

			imout.Pix[4*i] = uint8(k * float64(r))
			imout.Pix[4*i+1] = uint8(k * float64(g))
			imout.Pix[4*i+2] = uint8(k * float64(b))
			imout.Pix[4*i+3] = a
		}
		mmod = ImageProxy{Image: imout}
		mmod.AddMetadata("Equalised - 32-bit RGBA")
		fmt.Println("equalised RGBA")
		fmt.Println(imout.Bounds())
	case *image.NRGBA:
		// cheaty cheaty - take as v the max of r/g/b - store as a histogram, make the cdf and map the values from the array to these.
		im := ip.Image.(*image.NRGBA)
		w := im.Bounds().Dx()
		h := im.Bounds().Dy()
		N := w * h
		var r, g, b, a, v uint8
		mapper := make(map[int]int)
		for i := 0; i < N; i++ {
			r = im.Pix[4*i+0]
			g = im.Pix[4*i+1]
			b = im.Pix[4*i+2]
			v = max(r, g, b)
			mapper[int(v)]++
		}

		// our mapper is complete
		cdf := make([]float64, 256)
		cdf[0] = float64(mapper[0])
		for i := 1; i < 256; i++ {
			cdf[i] = float64(mapper[i]) + cdf[i-1]
		}
		// cdf complete
		// replace the mapper with the cdf function values
		for i := 0; i < 256; i++ {
			mapper[i] = int(cdf[i] * 255 / float64(N))
		}
		imout := image.NewNRGBA(im.Bounds())

		for i := 0; i < N; i++ {

			r = im.Pix[4*i+0]
			g = im.Pix[4*i+1]
			b = im.Pix[4*i+2]
			a = im.Pix[4*i+3]
			v1 := float64(max(r, g, b))
			// v1 is mapped to v2 by the equalisation
			v2 := float64(mapper[int(v1)])
			k := v2 / v1

			imout.Pix[4*i+0] = uint8(k * float64(r))
			imout.Pix[4*i+1] = uint8(k * float64(g))
			imout.Pix[4*i+2] = uint8(k * float64(b))
			imout.Pix[4*i+3] = a

		}

		mmod = ImageProxy{Image: imout}
		mmod.AddMetadata("Equalised - 32-bit NRGBA")
		fmt.Println("equalised NRGBA")
		fmt.Println(imout.Bounds())
	case *image.RGBA64:

		// cheaty cheaty - take as v the max of r/g/b - store as a histogram, make the cdf and map the values from the array to these.
		im := ip.Image.(*image.RGBA64)
		w := im.Bounds().Dx()
		h := im.Bounds().Dy()
		N := w * h
		var r, g, b, v uint16
		mapper := make(map[int]int)
		for i := 0; i < N; i++ {
			r = uint16(im.Pix[8*i+0])<<8 + uint16(im.Pix[8*i+1])
			g = uint16(im.Pix[8*i+2])<<8 + uint16(im.Pix[8*i+3])
			b = uint16(im.Pix[8*i+4])<<8 + uint16(im.Pix[8*i+5])
			v = max(r, g, b)
			mapper[int(v)]++
		}
		// our mapper is complete
		cdf := make([]float64, 65536)
		cdf[0] = float64(mapper[0])
		for i := 1; i < 65536; i++ {
			cdf[i] = float64(mapper[i]) + cdf[i-1]
		}
		// cdf complete
		// replace the mapper with the cdf function values
		for i := 0; i < 65536; i++ {
			mapper[i] = int(cdf[i] * 65535 / float64(N))
		}
		imout := image.NewRGBA64(im.Bounds())

		for i := 0; i < N; i++ {

			r = uint16(im.Pix[8*i+0])<<8 + uint16(im.Pix[8*i+1])
			g = uint16(im.Pix[8*i+2])<<8 + uint16(im.Pix[8*i+3])
			b = uint16(im.Pix[8*i+4])<<8 + uint16(im.Pix[8*i+5])

			v1 := float64(max(r, g, b))
			// v is mapped to v2 by the equalisation
			v2 := float64(mapper[int(v)])
			k := v2 / v1

			imout.Pix[8*i+0] = uint8(uint16(k*float64(r)) >> 8)
			imout.Pix[8*i+1] = uint8(uint16(k*float64(r)) & 255)
			imout.Pix[8*i+2] = uint8(uint16(k*float64(g)) >> 8)
			imout.Pix[8*i+3] = uint8(uint16(k*float64(g)) & 255)
			imout.Pix[8*i+4] = uint8(uint16(k*float64(b)) >> 8)
			imout.Pix[8*i+5] = uint8(uint16(k*float64(b)) & 255)

			imout.Pix[8*i+6] = im.Pix[8*i+6]
			imout.Pix[8*i+7] = im.Pix[8*i+7]
		}
		mmod = ImageProxy{Image: imout}
		mmod.AddMetadata("Equalised - 64-bit RGBA")
		fmt.Println("equalised RGBA64")
		fmt.Println(imout.Bounds())
	case *image.NRGBA64:

		// cheaty cheaty - take as v the max of r/g/b - store as a histogram, make the cdf and map the values from the array to these.
		im := ip.Image.(*image.NRGBA64)
		w := im.Bounds().Dx()
		h := im.Bounds().Dy()
		N := w * h
		var r, g, b, v uint16
		mapper := make(map[int]int)
		for i := 0; i < N; i++ {
			r = uint16(im.Pix[8*i+0])<<8 + uint16(im.Pix[8*i+1])
			g = uint16(im.Pix[8*i+2])<<8 + uint16(im.Pix[8*i+3])
			b = uint16(im.Pix[8*i+4])<<8 + uint16(im.Pix[8*i+5])
			v = max(r, g, b)
			mapper[int(v)]++
		}
		// our mapper is complete
		cdf := make([]float64, 65536)
		cdf[0] = float64(mapper[0])
		for i := 1; i < 65536; i++ {
			cdf[i] = float64(mapper[i]) + cdf[i-1]
		}
		// cdf complete
		// replace the mapper with the cdf function values
		for i := 0; i < 65536; i++ {
			mapper[i] = int(cdf[i] * 65535 / float64(N))
		}
		imout := image.NewNRGBA64(im.Bounds())

		for i := 0; i < N; i++ {

			r = uint16(im.Pix[8*i+0])<<8 + uint16(im.Pix[8*i+1])
			g = uint16(im.Pix[8*i+2])<<8 + uint16(im.Pix[8*i+3])
			b = uint16(im.Pix[8*i+4])<<8 + uint16(im.Pix[8*i+5])

			v1 := float64(max(r, g, b))
			// v is mapped to v2 by the equalisation
			v2 := float64(mapper[int(v)])
			k := v2 / v1

			imout.Pix[8*i+0] = uint8(uint16(k*float64(r)) >> 8)
			imout.Pix[8*i+1] = uint8(uint16(k*float64(r)) & 255)
			imout.Pix[8*i+2] = uint8(uint16(k*float64(g)) >> 8)
			imout.Pix[8*i+3] = uint8(uint16(k*float64(g)) & 255)
			imout.Pix[8*i+4] = uint8(uint16(k*float64(b)) >> 8)
			imout.Pix[8*i+5] = uint8(uint16(k*float64(b)) & 255)

			imout.Pix[8*i+6] = im.Pix[8*i+6]
			imout.Pix[8*i+7] = im.Pix[8*i+7]
		}
		mmod = ImageProxy{Image: imout}
		mmod.AddMetadata("Equalised - 64-bit NRGBA")
		fmt.Println("equalised NRGBA64")
		fmt.Println(imout.Bounds())

	default:
		fmt.Println("Unknown format")
	}

	return &mmod
}

func (ip ImageProxy) Rotate(angle float64) *ImageProxy {
	ip2 := new(ImageProxy)
	switch im := ip.Image.(type) {
	case *image.Gray:
		a := ip.LuminanceAsFloat().Rotate(angle).AsImage(8, fmt.Sprintf("Rotated by  %.2f radians", angle), false)
		return &a
	case *image.Gray16:
		a := ip.LuminanceAsFloat().Rotate(angle).AsImage(16, fmt.Sprintf("Rotated by  %.2f radians", angle), false)
		return &a
	case *image.RGBA:
		fmt.Println("32-bit RGBA")
		cs := make([]*Float2D, 4)
		// components in an array indexed by Component

		fmt.Println(R, G, B, A)

		cs[0] = ip.GetComponentAsFloat(R)
		fmt.Println(cs[0].Range())
		cs[0] = cs[0].Rotate(angle)
		fmt.Println(cs[0].Range())
		cs[1] = ip.GetComponentAsFloat(G)
		fmt.Println(cs[1].Range())
		cs[1] = cs[1].Rotate(angle)
		fmt.Println(cs[1].Range())
		cs[2] = ip.GetComponentAsFloat(B)
		fmt.Println(cs[2].Range())
		cs[2] = cs[2].Rotate(angle)
		fmt.Println(cs[2].Range())
		cs[3] = ip.GetComponentAsFloat(A)
		fmt.Println(cs[3].Range())
		cs[3] = cs[3].Rotate(angle)
		fmt.Println(cs[3].Range())

		// reassemble rgba
		h := len(*cs[0])
		w := len((*cs[0])[0])
		fmt.Println(w, h, "after rotation")
		// px := make([]uint8, 4*h*w)
		// // Pill the fixels
		// for j:=0;j<h;j++{
		// 	for i:=0;i<w;i++{
		// 		for k:=0;k<4;k++{
		// 			px[4*(j*w+i)+k]=uint8(cs[k][j][i])
		// 		}
		// 	}
		// }

		fmt.Println("RGBA")
		imout := image.NewRGBA(image.Rect(0, 0, w, h))

		// imout.Pix = px
		for x := 0; x < w; x++ {
			for y := 0; y < h; y++ {
				imout.SetRGBA(x, y, color.RGBA{uint8((*cs[0])[y][x]), uint8((*cs[1])[y][x]), uint8((*cs[2])[y][x]), 255})
			}
		}
		prox := ImageProxy{Image: imout}
		prox.metadata = ip.metadata
		prox.AddMetadata(fmt.Sprintf("Rotated %.2f radians RGBA", angle))
		return &prox

	case *image.NRGBA:
		fmt.Println("32-bit NRGBA")
		cs := make([]*Float2D, 4)
		// components in an array indexed by Component

		// fmt.Println(R, G, B, A)

		cs[0] = ip.GetComponentAsFloat(R)
		fmt.Println(cs[0].Range())
		cs[0] = cs[0].Rotate(angle)
		cs[1] = ip.GetComponentAsFloat(G)
		fmt.Println(cs[1].Range())
		cs[1] = cs[1].Rotate(angle)
		cs[2] = ip.GetComponentAsFloat(B)
		fmt.Println(cs[2].Range())
		cs[2] = cs[2].Rotate(angle)
		cs[3] = ip.GetComponentAsFloat(A)
		fmt.Println(cs[3].Range())
		cs[3] = cs[3].Rotate(angle)

		// reassemble rgba
		h := len(*cs[0])
		w := len((*cs[0])[0])
		fmt.Println(w, h, "after rotation")
		px := make([]uint8, 4*h*w)
		// Pill the fixels
		for j := 0; j < h; j++ {
			for i := 0; i < w; i++ {
				for k := 0; k < 4; k++ {
					px[4*(j*w+i)+k] = uint8((*cs[k])[j][i])
				}
			}
		}

		fmt.Println("RGBA")
		imout := image.NewRGBA(image.Rect(0, 0, w, h))

		imout.Pix = px

		prox := ImageProxy{Image: imout}
		prox.metadata = ip.metadata
		prox.AddMetadata(fmt.Sprintf("Rotated %.2f radians RGBA", angle))
		return &prox

	case *image.RGBA64, *image.NRGBA64:
		fmt.Println("64-bit")
		cs := make([]*Float2D, 4)
		// components in an array indexed by Component
		cs[0] = ip.GetComponentAsFloat(R).Rotate(angle)
		cs[1] = ip.GetComponentAsFloat(G).Rotate(angle)
		cs[2] = ip.GetComponentAsFloat(B).Rotate(angle)
		cs[3] = ip.GetComponentAsFloat(A).Rotate(angle)
		// reassemble rgba
		h := len(*cs[0])
		w := len((*cs[0])[0])
		px := make([]uint8, 8*h*w)
		// Pill the fixels
		for j := 0; j < h; j++ {
			for i := 0; i < w; i++ {
				for c := 0; c < 4; c++ {
					k := uint16((*cs[c])[j][i])
					px[8*(j*w+i)+2*c] = uint8(k >> 8)
					px[8*(j*w+i)+2*c+1] = uint8(k & 255)
				}
			}
		}
		switch im.(type) {
		case *image.RGBA64:
			imout := image.NewRGBA64(image.Rect(0, 0, w, h))
			imout.Pix = px
			prox := ImageProxy{Image: imout}
			prox.metadata = ip.metadata
			prox.AddMetadata(fmt.Sprintf("Rotated %.2f radians", angle))
			return &prox
		case *image.NRGBA64:
			imout := image.NewNRGBA64(image.Rect(0, 0, w, h))
			imout.Pix = px
			prox := ImageProxy{Image: imout}
			prox.metadata = ip.metadata
			prox.AddMetadata(fmt.Sprintf("Rotated %.2f radians", angle))
			return &prox
		}
	}

	return ip2
}

// modified rotator using existing Float2D functionality
func (ip ImageProxy) Rotate2(angle float64) *ImageProxy {

	newim := imaging.Rotate(ip, angle*180/math.Pi, color.Transparent)

	op := new(ImageProxy)
	op.Image = newim
	op.AddMetadata(fmt.Sprintf("Gray Rotated by  %.2f radians", angle))
	return op
}

// Fourier-domain sub-pixel shift
//
// Returned image is the same size as the original.
//
//	x,y - amount by which to shift, in pixels.
func (ip ImageProxy) TranslateSubpixel(x, y float64) *ImageProxy {
	// h := ip.Bounds().Dy()
	w := ip.Bounds().Dx()

	switch im := ip.Image.(type) {
	case *image.Gray:
		log.Println("GRAY")
		a := ListToGrid(im.Pix, w).AsFloat().TranslateSubpixel(x, y).AsImage(8, fmt.Sprintf("TX: %.2f,%.2f", x, y), false)
		return &a

	case *image.Gray16:
		px := Short8As16(im.Pix)
		a := ListToGrid(px, w).AsFloat().TranslateSubpixel(x, y).AsImage(8, fmt.Sprintf("TX: %.2f,%.2f", x, y), false)
		return &a

	case *image.RGBA, *image.NRGBA:
		ipout := ImageProxy{}
		log.Println("RGBA 32 Bit")
		rotatedComponents := make([]Float2D, 4)
		wg := &sync.WaitGroup{}
		wg.Add(3)
		for c := range []int{0, 1, 2} {
			go func(c int, wg *sync.WaitGroup) {
				defer wg.Done()
				rotatedComponents[c] = ip.GetComponentAsFloat(Component(c)).TranslateSubpixel(x, y)
			}(c, wg)
		}
		wg.Wait()

		ipout.CreateFromRGB(rotatedComponents[:3], true)

		// fmt.Println(a.Bounds())
		return &ipout

	case *image.RGBA64, *image.NRGBA64:

		rotatedComponents := make([]Float2D, 4)
		for c := range []int{0, 1, 2, 3} {
			rotatedComponents[c] = ip.GetComponentAsFloat(Component(c)).TranslateSubpixel(x, y)
		}
		// all 4 components have been rotated - recombine as RGBA
		a := ImageProxy{}
		a.LoadFromFloatsAsRGBAComponents(rotatedComponents, 16, false)
		return &a
	default:
		return &ip
	}
}

func (ip *ImageProxy) MakePyramid(minsize int) ([]*ImageProxy, error) {
	h := ip.Bounds().Dx()
	if ip.Bounds().Dy() < h {
		h = ip.Bounds().Dy()
	}
	if h < minsize {
		return nil, errors.New("image is already smaller than the required minimum")
	}

	var layers []*ImageProxy = []*ImageProxy{ip}

	for h > minsize {
		lastlayer := layers[len(layers)-1]
		layers = append(layers, lastlayer.halveImage())
		h /= 2
	}
	return layers, nil
}
func (ip *ImageProxy) halveImage() *ImageProxy {
	b := ip.Bounds()
	op := new(ImageProxy)
	op.Image = imaging.Resize(ip, b.Dx()/2, b.Dy()/2, imaging.Gaussian)
	op.AddMetadata("Gaussian resample + reduction")
	return op
}
