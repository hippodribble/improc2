package proxy


import (
	_ "image/jpeg"
	"image/png"
	_ "image/png"
	"log"
	"math"
	"os"
	"testing"
)

// func TestLoadProxyImage(t *testing.T) {
// 	ip := new(ImageProxy)
// 	err := ip.LoadFromFile("testdata/Como Blue.png")
// 	if err != nil {
// 		t.Error(err)
// 	}
// 	if ip.Image.Bounds().Dx() <= 0 {
// 		t.Error("Loaded image is of invalid size")
// 	}
// }

// func TestLoadProxyImageFromFloat2D(t *testing.T) {
// 	ip := new(ImageProxy)
// 	floatIn := *NewFloat2D(10, 10)
// 	floatIn[2][1] = 300
// 	floatIn[2][2] = -100
// 	err := ip.LoadFromFloats(floatIn, 8, false)
// 	if err != nil {
// 		t.Error(err)
// 	}
// 	if ip.Image.Bounds().Dx() <= 0 {
// 		t.Error("Loaded image is of invalid size")
// 	}
// 	// for j:=0;j<10;j++{
// 	// 	s:=""
// 	// 	for i:=0;i<10;i++{
// 	// 		s+=fmt.Sprintf("%05d",ip.At(i,j))
// 	// 	}
// 	// 	log.Println(s)
// 	// }

// }

// func TestLoadProxyImageFromFloatsAsRGBAComponents(t *testing.T) {
// 	ip := new(ImageProxy)
// 	r := *NewFloat2D(10, 10)
// 	g := *NewFloat2D(10, 10)
// 	b := *NewFloat2D(10, 10)
// 	a := *NewFloat2D(10, 10)
// 	for j := 0; j < 10; j++ {
// 		for i := 0; i < 10; i++ {
// 			r[j][i] = float64(i * j)
// 			g[j][i] = float64((10 - i) * j)
// 			b[j][i] = float64(i * (10 - j))
// 			a[j][i] = float64(255)
// 		}
// 	}
// 	err := ip.LoadFromFloatsAsRGBAComponents([]Float2D{r, g, b, a}, 16, false)
// 	if err != nil {
// 		t.Error(err)
// 	}
// 	if ip.Image.Bounds().Dx() <= 0 {
// 		t.Error("Loaded image is of invalid size")
// 	}
// }

// func TestProxyImageGetComponent(t *testing.T) {

// 	// err := ip.LoadFromFile("testdata/Como Blue.png")
// 	// if err != nil {
// 	// 	t.Error(err)
// 	// }
// 	ip := new(ImageProxy)
// 	ip.Format = "memory"
// 	ip.Config = image.Config{ColorModel: color.NRGBAModel, Width: 5, Height: 5}
// 	ip.Image = image.NewRGBA(image.Rect(0, 0, 5, 5))
// 	for j := 0; j < ip.Config.Height; j++ {
// 		for i := 0; i < ip.Config.Width; i++ {
// 			ip.Image.(*image.RGBA).SetRGBA(i, j, color.RGBA{R: 1, G: 2, B: 4, A: 8})
// 		}
// 	}
// 	temp := *ip.GetComponent(R)
// 	r := temp.Image.(*image.Gray)
// 	temp = *ip.GetComponent(G)
// 	g := temp.Image.(*image.Gray)
// 	temp = *ip.GetComponent(B)
// 	b := temp.Image.(*image.Gray)
// 	temp = *ip.GetComponent(A)
// 	a := temp.Image.(*image.Gray)

// 	if r.Bounds() != ip.Bounds() {
// 		t.Error("lengths of component and image pix arrays were not equal")
// 	}
// 	for j := 0; j < ip.Bounds().Dy(); j++ {
// 		for i := 0; i < ip.Bounds().Dx(); i++ {
// 			if r.Pix[j*ip.Config.Width+i] != uint8(1) {
// 				t.Error("component R had unexpected value: " + fmt.Sprintf("%d,%d %v", i, j, r.Pix[j*ip.Config.Width+i]))
// 			}
// 		}
// 	}

// 	if g.Bounds() != ip.Bounds() {
// 		t.Error("lengths of component and image pix arrays were not equal")
// 	}
// 	for j := 0; j < ip.Bounds().Dy(); j++ {
// 		for i := 0; i < ip.Bounds().Dx(); i++ {
// 			if g.Pix[j*ip.Config.Width+i] != uint8(2) {
// 				t.Error("component G had unexpected value: " + fmt.Sprintf("%d,%d %v", i, j, g.Pix[j*ip.Config.Width+i]))
// 			}
// 		}
// 	}

// 	if b.Bounds() != ip.Bounds() {
// 		t.Error("lengths of component and image pix arrays were not equal")
// 	}
// 	for j := 0; j < ip.Bounds().Dy(); j++ {
// 		for i := 0; i < ip.Bounds().Dx(); i++ {
// 			if b.Pix[j*ip.Config.Width+i] != uint8(4) {
// 				t.Error("component B had unexpected value: " + fmt.Sprintf("%d,%d %v", i, j, b.Pix[j*ip.Config.Width+i]))
// 			}
// 		}
// 	}

// 	if a.Bounds() != ip.Bounds() {
// 		t.Error("lengths of component and image pix arrays were not equal")
// 	}
// 	for j := 0; j < ip.Bounds().Dy(); j++ {
// 		for i := 0; i < ip.Bounds().Dx(); i++ {
// 			if a.Pix[j*ip.Config.Width+i] != uint8(8) {
// 				t.Error("component A had unexpected value: " + fmt.Sprintf("%d,%d %v", i, j, a.Pix[j*ip.Config.Width+i]))
// 			}
// 		}
// 	}

// }

// func TestImageProxyCopy(t *testing.T) {
// 	ip := new(ImageProxy)
// 	err := ip.LoadFromFile("testdata/Como Blue.png")
// 	if err != nil {
// 		t.Error("File not opened")
// 	}
// 	t.Log("Start Copy")
// 	ip2 := ip.Copy()
// 	t.Log("Copy finished")

// 	f, err := os.Create("testdata/copy.png")
// 	if err != nil {
// 		t.Error(err)
// 	}
// 	defer f.Close()
// 	png.Encode(f, ip2)
// }

// func TestImageProxyCircularise(t *testing.T) {
// 	ip := new(ImageProxy)
// 	err := ip.LoadFromFile("testdata/Como Blue.png")
// 	if err != nil {
// 		t.Error("File not opened")
// 	}

// 	t.Log("Start Circularise")
// 	c := ip.Circularise()
// 	t.Log("Circle complete")

// 	g, err := os.Create("testdata/circle.png")
// 	if err != nil {
// 		t.Error(err)
// 	}
// 	defer g.Close()
// 	png.Encode(g, c)

// }

// func TestWriteImage(t *testing.T){
// 	t.Logf("Loading at %s\n",time.Now())
// 	ip := new(ImageProxy)
// 	err := ip.LoadFromFile("testdata/Como Blue.png")
// 	if err != nil {
// 		t.Error("File not opened")
// 	}
// 	t.Logf("Loaded at %s\n",time.Now())
// 	t.Logf("Writing at %s\n",time.Now())

// 	f,err:=os.Create("testdata/writespeedcheck.png")
// 	if err!=nil{
// 		t.Error("Failed to open file for writing")
// 	}
// 	defer f.Close()
// 	err=png.Encode(f,ip)
// 	if err!=nil{
// 		t.Error("Failed to encode image")
// 	}
// 	t.Logf("Complete at %s\n",time.Now())

// }

// func TestHistogramEqualise(t *testing.T) {
// 	ip := new(ImageProxy)
// 	err := ip.LoadFromFile("testdata/Como Blue.png")
// 	if err != nil {
// 		t.Error("File not opened")
// 	}
// 	ip2 := ip.HistogramEqualise()
// 	g, err := os.Create("testdata/histogramequalised.png")
// 	if err != nil {
// 		t.Error(err)
// 	}
// 	defer g.Close()
// 	png.Encode(g, ip2)
// }

func TestImageProxyRotate(t *testing.T) {
	ip := new(ImageProxy)
	err := ip.LoadFromFile("testdata/IMG_20240823_144813_EXP3.jpg")
	// err := ip.LoadFromFile("testdata/indicator.png")
	if err != nil {
		t.Error("File not opened")
	}

	log.Println("after loading", ip.Type(), ip.ColorModel())
	ip3 := ip.Rotate2(math.Pi / 4)
	// ip3 := ip.Rotate2(0)
	// ip3:=ip.Copy()
	log.Println(ip3.AllMetadata())
	h, err := os.Create("testdata/rotated2.png")
	if err != nil {
		t.Error(err)
	}
	defer h.Close()
	png.Encode(h, ip3)
}
