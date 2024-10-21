package proxy



func equalise8bit(pix *[]uint8) *[]uint8 {

	pxout := make([]uint8, len(*pix))
	copy(pxout, *pix)
	mapper := make(map[int]int)
	// warm up the map
	for i := 0; i < 256; i++ {
		mapper[i] = 0
	}
	N := len(*pix)

	// make histogram
	for i := 0; i < N; i++ {
		I := int((*pix)[i])
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

	for i := 0; i < len(*pix); i++ {
		index := int((*pix)[i])
		pxout[i] = uint8(mapper[index])
	}
	return &pxout
}
func equalise16bit(pix *[]uint8) *[]uint8 {

	pxout := make([]uint8, len(*pix))
	mapper := make(map[int]int)
	// warm up the map
	for i := 0; i < 65536; i++ {
		mapper[i] = 0
	}
	N := len(*pix) / 2

	// make histogram
	for i := 0; i < N; i++ {
		I := int((*pix)[2*i])*256 + int((*pix)[2*i+1])
		// I := int((*pix)[2*i]*256+(*pix)[2*i+1])
		mapper[I]++
	}
	
	cdf := make([]float64, 65536)
	cdf[0] = float64(mapper[0])
	// var j uint8
	for i := 1; i < 65536; i++ {
		cdf[i] = cdf[i-1] + float64(mapper[i])
	}

	for i := 0; i < 65536; i++ {
		mapper[i] = int(cdf[i] * 65536 / float64(N))
	}

	for i := 0; i < len((*pix))/2; i++ {
		index := int((*pix)[2*i])*256 + int((*pix)[2*i+1])
		pxout[2*i] = uint8((mapper[index]) >> 8)
		pxout[2*i+1] = uint8((mapper[index] & 255))
	}
	return &pxout
}

func equalise16(v *[]uint16)*[]uint16{
	N:=len(*v)
	vout := make([]uint16, N)
	mapper := make(map[int]int)
	// warm up the map
	for i := 0; i < 65536; i++ {
		mapper[i] = 0
	}
	// make histogram
	for i := 0; i < N; i++ {
		I := (*v)[i]
		// I := int((*pix)[2*i]*256+(*pix)[2*i+1])
		mapper[int(I)]++
	}
	
	cdf := make([]float64, 65536)
	cdf[0] = float64(mapper[0])
	// var j uint8
	for i := 1; i < 65536; i++ {
		cdf[i] = cdf[i-1] + float64(mapper[i])
	}

	for i := 0; i < 65536; i++ {
		mapper[i] = int(cdf[i] * 65536 / float64(N))
	}
	
	for i := 0; i < N; i++ {
		index := int((*v)[i])
		vout[i] = uint16(mapper[index])
	}
	return &vout

}

func convert8to16(in *[]uint8) *[]uint16{
	out := make([]uint16, len(*in)/2)
	for i := 0; i < len(*in)/2; i++ {
		out[i] = uint16((*in)[2*i])<<8 | uint16((*in)[2*i+1])
	}
	return &out

}
func convert16to8(in *[]uint16) *[]uint8{
	out := make([]uint8, len(*in)*2)
	for i := 0; i < len(*in); i++ {
		out[2*i] = uint8((*in)[i] >> 8)
		out[2*i+1] = uint8((*in)[i] & 255)
	}
	return &out
}