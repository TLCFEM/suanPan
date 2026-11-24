// Minimal declaration of the oneMKL RCI FGMRES entry points used by the
// standalone tests.  The real headers ship with oneMKL; declaring only what the
// wrapper needs keeps the test applications self-contained.
#ifndef TEST_FGMRES_MKL_RCI_H
#define TEST_FGMRES_MKL_RCI_H

using MKL_INT = int;

extern "C" {
void dfgmres_init(const MKL_INT* n, const double* x, const double* b, MKL_INT* RCI_request, MKL_INT* ipar, double* dpar, double* tmp);
void dfgmres_check(const MKL_INT* n, const double* x, const double* b, MKL_INT* RCI_request, MKL_INT* ipar, double* dpar, double* tmp);
void dfgmres(const MKL_INT* n, double* x, double* b, MKL_INT* RCI_request, MKL_INT* ipar, double* dpar, double* tmp);
void dfgmres_get(const MKL_INT* n, double* x, double* b, MKL_INT* RCI_request, const MKL_INT* ipar, const double* dpar, double* tmp, MKL_INT* itercount);
}

#endif
