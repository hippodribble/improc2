package main

import (
	"errors"
	"fmt"
	"image"
	"image/color"
	"log"
	"os"
	"time"

	"fyne.io/fyne/v2"
	"fyne.io/fyne/v2/app"
	"fyne.io/fyne/v2/canvas"
	"fyne.io/fyne/v2/container"
	"fyne.io/fyne/v2/widget"
	"github.com/disintegration/imaging"
	"github.com/hippodribble/fynewidgets"
	"github.com/hippodribble/improc2/proxy"
)

var stack *fyne.Container
var win fyne.Window
var status *widget.Label
var levelslider *widget.Slider
var ci *canvas.Image
var infochan chan interface{}

// var zc ZoomCanvas
var pyramid []*proxy.ImageProxy

func main() {
	log.Println("Application start")
	win = simpleGUI()
	// win.SetFullScreen(true)
	win.Resize(fyne.NewSize(500, 500))
	win.ShowAndRun()
}

func simpleGUI() fyne.Window {
	app := app.NewWithID("com.github.hippodribble.improc2")
	w := app.NewWindow("Image Fusion")
	stack = container.NewStack(widget.NewLabel("Welcome"))
	top := container.NewHBox()
	levelslider = widget.NewSlider(0, 0)
	levelslider.Step = 1
	levelslider.OnChanged = setimageindex
	r := canvas.NewRectangle(color.Transparent)
	r.SetMinSize(fyne.NewSize(300, 5))
	a := container.NewStack(r, levelslider)
	top.Add(a)
	top.Add(widget.NewSeparator())
	top.Add(widget.NewButton("Single", makeSingle))
	top.Add(widget.NewButton("Info", func() {
		go func() {
			w := fyne.CurrentApp().Driver().AllWindows()[0]
			w.Show()
			w.SetFullScreen(true)
			sz := w.Canvas().Size()
			w.SetFullScreen(false)
			log.Println(sz)
		}()
	}))

	infochan = make(chan interface{})
	bottom := fynewidgets.NewStatusProgress(infochan)

	w.SetContent(container.NewBorder(top, bottom, nil, nil, stack))
	return w
}

func setimageindex(f float64) {
	if pyramid == nil {
		return
	}
	if ci == nil {
		return
	}
	i := int(f)
	if i < len(pyramid) {
		status.SetText(fmt.Sprintf("zoom level %v selected", i))
		ci.Image = pyramid[i]
		ci.Refresh()

		// zc.Image.Image = pyramid[i]
		// zc.Refresh()
	}
}

func makeSingle() {
	name := "proxy/testdata/backup/200123_Maximum resolution.jpg"

	f, err := os.Open(name)

	if err != nil {
		log.Println("error reading file " + err.Error())
		return
	}

	ip, _, err := image.Decode(f)
	if err != nil {
		log.Println("could not decode image " + err.Error())
	}

	// ww, err := NewImageWidget(ip, 100)
	ww, err := fynewidgets.NewPanZoomCanvasFromImage(ip, image.Pt(100, 100), infochan, name)
	if err != nil {
		log.Println("cannot make pyramid from images")
		return
	}
	infochan <- fmt.Sprintf("File %s (level %d)", name, ww.Datum().Pyramid.Level())
	stack.RemoveAll()
	stack.Add(ww)
	stack.Refresh()
}

type ImageWidget struct {
	widget.BaseWidget
	Image      canvas.Image
	Pyramid    []*image.Image
	Index      int
	pause      bool
	updateRate int `default:"1000"`
}

func NewImageWidget(img image.Image, minsize int) (*ImageWidget, error) {
	pyr, err := makePyramid(&img, minsize)
	if err != nil {
		return nil, errors.Join(err)
	}
	index := len(pyr) - 1
	ci := *canvas.NewImageFromImage(*pyr[index])
	ci.FillMode = canvas.ImageFillContain
	w := &ImageWidget{
		Image:   ci,
		Pyramid: pyr,
		Index:   index,
	}
	w.ExtendBaseWidget(w)

	var xratio, yratio, ratio float32

	go func() {
		for {
			time.Sleep(time.Millisecond * time.Duration(w.updateRate))
			if w.pause {
				continue
			}

			size := w.Size()
			xratio = size.Width / float32((*w.Pyramid[w.Index]).Bounds().Dx())
			yratio = size.Height / float32((*w.Pyramid[w.Index]).Bounds().Dy())
			ratio = min(xratio, yratio)

			if ratio < 1.1 && ratio > .5 {
				continue
			}
			if ratio > 1.1 {
				w.higherResolution()
			} else {
				w.lowerResolution()
			}
		}
	}()
	return w, nil
}

func (item *ImageWidget) CreateRenderer() fyne.WidgetRenderer {
	c := container.NewStack(&item.Image)
	ren := widget.NewSimpleRenderer(c)
	return ren
}

func (item *ImageWidget) SetUpdateRate(milliseconds int) error {
	if milliseconds < 1 || milliseconds > 10000 {
		return errors.New("update rate should be between 1 and 10000")
	}
	item.updateRate = milliseconds
	return nil
}

func (m *ImageWidget) higherResolution() {
	if m.Index == 0 {
		return
	}
	m.Index -= 1
	m.changeResolution()
}

func (m *ImageWidget) lowerResolution() {
	if m.Index == len(m.Pyramid)-1 {
		return
	}
	m.Index += 1
	m.changeResolution()
}

func (m *ImageWidget) changeResolution() {
	ci := canvas.NewImageFromImage(*m.Pyramid[m.Index])
	ci.FillMode = canvas.ImageFillContain
	m.Image = *ci
	m.Refresh()
}

func makePyramid(img *image.Image, minsize int) ([]*image.Image, error) {
	h := (*img).Bounds().Dx()
	if (*img).Bounds().Dy() < h {
		h = (*img).Bounds().Dy()
	}
	if h < minsize {
		return nil, errors.New("image is already smaller than the required minimum")
	}

	var layers []*image.Image = []*image.Image{img}

	for h > minsize {
		lastlayer := layers[len(layers)-1]

		b := (*lastlayer).Bounds()
		newlayer := image.Image(imaging.Resize(*lastlayer, b.Dx()/2, b.Dy()/2, imaging.Gaussian))
		layers = append(layers, &newlayer)
		h /= 2
	}
	return layers, nil
}
