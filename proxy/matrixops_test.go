package proxy

import (
	"log"
	"testing"

	"gonum.org/v1/gonum/mat"
)

func TestDeterminant(t *testing.T) {

	m := NewImageMatrix(mat.NewDense(2, 2, []float64{1, 2, 3, 4}))
	detM := m.Determinant()
	log.Println(detM)

}
