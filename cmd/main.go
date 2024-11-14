package main

import (
	"errors"
	"fmt"
	"image"
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
	// stack = container.NewStack(widget.NewLabel("Welcome"))
	stack = container.NewGridWithColumns(1, widget.NewLabel("Welcome"))
	
	top := container.NewHBox()

	// r := canvas.NewRectangle(color.Transparent)
	// r.SetMinSize(fyne.NewSize(300, 5))

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
	w.SetContent(container.NewBorder(top, bottom, nil, nil, container.NewVScroll(stack)))
	return w
}

func makeSingle() {
	name := "/Users/glenn/Dropbox/golang/imagewidgets/moon.png"

	f, err := os.Open(name)

	if err != nil {
		log.Println("error reading file " + err.Error())
		return
	}

	ip, _, err := image.Decode(f)
	if err != nil {
		log.Println("could not decode image " + err.Error())
	}

	ww, err := fynewidgets.NewPanZoomCanvasFromImage(ip, image.Pt(100, 100), infochan, name)
	if err != nil {
		log.Println("cannot make pyramid from images")
		return
	}

	infochan <- fmt.Sprintf("File %s (level %d)", name, ww.Datum().Pyramid.Level())
	stack.RemoveAll()
	stack.Add(ww)
	stack.Refresh()

	improx:=&proxy.ImageProxy{Image: ip}
	pf:=improx.LuminanceAsFloat()
	log.Println("luminance")

	filtered:=pf.BilateralFilter(2,2,2)
	prox2:=&proxy.ImageProxy{}
	prox2.LoadFromFloats(filtered,8,true)
	w2,_:=fynewidgets.NewPanZoomCanvasFromImage(prox2.Image, image.Pt(100, 100), infochan, name+" bilateral")
	stack.Add(w2)
	stack.Refresh()

	// filtered2:=pf.BilateralFilter(5,5,5)
	// prox3:=&proxy.ImageProxy{}
	// prox3.LoadFromFloats(filtered2,8,true)
	// w3,_:=fynewidgets.NewPanZoomCanvasFromImage(prox3.Image, image.Pt(100, 100), infochan, name+" bilateral 2")
	// stack.Add(w3)
	// stack.Refresh()
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
