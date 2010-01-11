;;; -*- Mode: Lisp; Package: BORDEAUX-FFT -*-

;;;  (c) copyright 2004-2009 by
;;;           Robert Strandh <strandh@labri.fr>
;;;           Sylvain Marchand <sm@labri.fr>
;;;           Martin Raspaud <zorgzorg2@gmail.com>
;;;  (c) copyright 2008-2009 by
;;;           Andy Hefner <ahefner@gmail.com>
;;;  (c) copyright 2009 by
;;;           Paul Khuong <pvk@pvk.ca>
;;;
;;; This program is free software; you can redistribute it and/or
;;; modify it under the terms of the GNU General Public License
;;; as published by the Free Software Foundation; either version 2
;;; of the License, or (at your option) any later version.
;;;
;;; This program is distributed in the hope that it will be useful,
;;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;;; GNU General Public License for more details.
;;;
;;; You should have received a copy of the GNU General Public License
;;; along with this program; if not, write to the
;;; Free Software Foundation, Inc., 59 Temple Place - Suite 330,
;;; Boston, MA 02111-1307, USA.

(defpackage :bordeaux-fft 
  (:use :common-lisp)
  (:export :fft :sfft :fft!
           :ifft :sifft :ifft!
           :complex-sample :complex-sample-array
           :rectangular 
           :hann
           :blackman
           :triangle
           :bartlett
           :gaussian
           :gaussian*bartlett^x
           :blackman-harris
           :window-vector
           :extract-window
           :extract-window-into
           :extract-centered-window
           :extract-centered-window-into
           :windowed-fft
           :*fft-instance*
           :*ifft-instance*))

(in-package :bordeaux-fft)

(deftype complex-sample () `(complex double-float))
(deftype complex-sample-array ()  `(simple-array complex-sample (*)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; FFT implementation

(defclass fft-instance ()
  ((size   :initarg :size :reader size)
   (inter  :initarg :inter)
   (coeffs :initarg :coeffs)))

(defun make-fourier-instance (size direction) 
  ;; direction 1 direct, -1 inverse
  ;; check here that size is a power of two
  (let ((inter (make-array size :element-type 'complex-sample))
        (coeffs (make-array size :element-type 'complex-sample))
        (coeff #c(1d0 0d0))
        (factor (exp (/ (* direction -2 pi #c(0 1)) size))))
    (loop for i of-type fixnum from (ash size -1) below size
          do (setf (aref coeffs i) coeff)
          do (setf coeff (* coeff factor)))
    (do ((src (ash size -1) (ash src -1))
         (src-size (ash size -1) (ash src-size -1))
         (dst (ash size -2) (ash dst -1)))
        ((< src-size 1))
      (loop for s from src below (+ src src-size) by 2
            for d upfrom dst
            do (setf (aref coeffs d) (aref coeffs s))))
    (make-instance 'fft-instance
                   :size size
                   :inter inter
                   :coeffs coeffs)))

(defun make-fft-instance (size)
  (make-fourier-instance size 1))

(defun fft-common (instance source dest)
  (let ((inter (slot-value instance 'inter))
        (coeffs (slot-value instance 'coeffs))
        (size (slot-value instance 'size)))
    (declare (type complex-sample-array
		   source dest inter coeffs))
    (declare (type (and unsigned-byte fixnum) size))
    (assert (= size (length source)))
    (assert (= size (length dest)))
    (labels ((aux (inter dest n starts steps startd)
               (declare (optimize (speed 3) (safety 0))
                        (type complex-sample-array inter dest)
                        (type (and fixnum unsigned-byte) starts steps startd))
               (if (= n 2)
                   (let ((a (aref source starts))
                         (b (aref source (+ starts steps))))
                     (declare (type (complex double-float) a b))
                     (setf (aref dest startd) (+ a b)
                           (aref dest (1+ startd)) (- a b))
                     nil)
                   (let ((2*steps (ash steps 1))
                         (n/2 (ash n -1)))
                     (declare (type fixnum 2*steps n/2))
                     (aux dest inter
                          n/2 starts 2*steps startd)
                     (aux dest inter n/2 
			  (+ starts steps) 
			  2*steps 
			  (+ startd n/2))
                     (loop for i of-type fixnum 
			   from (the fixnum (+ startd n/2)) by 2
                           for c of-type fixnum from n/2 by 2
                           for dummy of-type fixnum from 0 below (truncate n/2 2)
                           do (let ((i0 (aref inter i))
                                    (i1 (aref inter (1+ i)))
                                    (c0 (aref coeffs c))
                                    (c1 (aref coeffs (1+ c))))
                                (setf (aref inter i)      (* i0 c0)
                                      (aref inter (1+ i)) (* i1 c1))))
                     (loop for i of-type fixnum from startd by 2
                           for j of-type fixnum 
			   from (the fixnum (+ startd n/2)) by 2
                           for dummy of-type fixnum from 0 below (truncate n/2 2)
                           do (let ((a0 (aref inter i))
                                    (a1 (aref inter (1+ i)))
                                    (b0 (aref inter j))
                                    (b1 (aref inter (1+ j))))
                                (setf (aref dest i)      (+ a0 b0)
                                      (aref dest (1+ j)) (- a1 b1)
                                      (aref dest (1+ i)) (+ a1 b1)
                                      (aref dest j)      (- a0 b0))))
                     nil))))
      (aux inter dest size 0 1 0)
      dest)))

(defvar *fft-instance* (make-fft-instance 1024))

(defun fft (source)
  "Returns the Fourier transform of source, allocating a new array
   for the result."
  (unless (and *fft-instance*
               (= (size *fft-instance*) (length source)))
    (setf *fft-instance* (make-fft-instance (length source))))
  (let ((dest (make-array (length source)
                          :element-type '(complex double-float))))
    (fft-common *fft-instance* source dest)))
  
(defun fft! (source dest)
  "Destructive version of fft, since it fills dest."
  (unless (and *fft-instance*
               (= (size *fft-instance*) (length source)))
    (setf *fft-instance* (make-fft-instance (length source))))
  (fft-common *fft-instance* source dest))

(defun power-of-two (x)
  (and (> x 0)
       (zerop (logand x (1- x)))
       (1- (integer-length x))))

(defun sfft (source &optional len)
  "This is the generic fft function. Stands for stupid fft. Can take
   any kind of array as input."
  (let* ((l (or len (length source)))
         (power (power-of-two l))
	 (lp (expt 2 (or power (error "Input size is not a power of two"))))
	 (new-source (make-array lp :element-type 'complex-sample)))
    (loop for i fixnum from 0 below l do
      (setf (aref new-source i)
            (complex (float (realpart (aref source i)) 0.0d0)
                     (float (imagpart (aref source i)) 0.0d0))))
    (fft new-source)))

(defmacro with-fft-instance (instance size &body body)
  `(let ((,instance (make-ftt-instance ,size)))
    ,@body))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Inverse Fast Fourier Transform

(defun make-ifft-instance (size)
  (make-fourier-instance size -1))

(defun ifft-common (instance source dest)
  (declare (type complex-sample-array dest))
  (fft-common instance source dest)
  (let ((1/length (/ (float (length source) 0.0d0))))
    (map-into dest (lambda (x) (* x 1/length )) dest)))

(defvar *ifft-instance* (make-ifft-instance 1024))

(defun ifft (source)
  "Returns the inverse Fourier transform of source, allocating a new
   array for the result."
  (unless (and *ifft-instance*
               (= (size *ifft-instance*) (length source)))
    (setf *ifft-instance* (make-ifft-instance (length source))))
  (let ((dest (make-array (length source)
                          :element-type '(complex double-float))))
    (ifft-common *ifft-instance* source dest)))

(defun ifft! (source dest)
  "Destructive version of ifft, since it fills dest."
  (unless (and *ifft-instance*
               (= (size *ifft-instance*) (length source)))
    (setf *ifft-instance* (make-ifft-instance (length source))))
  (ifft-common *ifft-instance* source dest))

(defun sifft (source)
  "This is the generic fft function. Stands for stupid fft. Can take
   any kind of array as input."
  (let* ((l (length source))
	 (lp (expt 2 (power-of-two l)))
	 (new-source (if (typep source 'complex-sample-array)
			 source
			 (map 'complex-sample-array 
			      #'(lambda (x)
				  (coerce x 'complex-sample))
			      source)))
	 (modified-source
	  (adjust-array new-source lp :initial-element #c(0.0d0 0.0d0)
			:element-type 'complex-sample)))
    (ifft modified-source)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Signal windows

(defun rectangular (i n)
  (declare (ignore i n))
  1.0f0)

(defun hann (i n) (* 0.5 (- 1.0 (cos (/ (* 2 pi i) (1- n))))))

(defun blackman* (alpha i n)
  (let ((a0 (/ (- 1 alpha) 2))
        (a1 0.5)
        (a2 (/ alpha 2)))
    (+ a0 
       (- (* a1 (cos (/ (* 2 pi i) (1- n))))) 
       (* a2 (cos (/ (* 4 pi i) (1- n)))))))

(defun blackman (i n) (blackman* 0.16 i n))

(defun triangle (i n) 
  (* (/ 2 n) (- (* n 0.5) (abs (- i (* 0.5 (1- n)))))))

(defun bartlett (i n) 
  (* (/ 2 (1- n)) (- (* (1- n) 0.5) (abs (- i (* 0.5 (1- n)))))))

(defun gauss* (sigma i n)
  (let (([n-1]/2 (* 0.5 (1- n))))
    (exp (* -0.5 (expt (/ (- i [n-1]/2) (* sigma [n-1]/2)) 2)))))

(let ((cache (make-hash-table)))
  (defun gaussian (sigma) 
    (or (gethash sigma cache)
        (setf (gethash sigma cache)
              (lambda (i n) (gauss* sigma i n))))))

(let ((cache (make-hash-table :test 'equal)))
  (defun gaussian*bartlett^x (sigma triangle-exponent)
    (or (gethash (list sigma triangle-exponent) cache)
        (setf (gethash (list sigma triangle-exponent) cache)
              (lambda (i n)
                (* (realpart (expt (bartlett i n) triangle-exponent))
                   (gauss* sigma i n)))))))

(defun cosine-series (i n a0 a1 a2 a3)
  (flet ((f (scale x) (* scale (cos (/ (* x pi i) (1- n))))))
    (+ a0 (- (f a1 2)) (f a2 4) (- (f a3 6)))))

(defun blackman-harris (i n)
  (cosine-series i n 0.35875f0 0.48829f0 0.14128f0 0.01168f0))

(let ((cache (make-hash-table :test 'equalp)))
  (defun window-vector (function n)
    (or (gethash (list function n) cache)
        (setf (gethash (list function n) cache)
              (let ((v (make-sequence '(simple-array double-float (*)) n)))
                (dotimes (i n v) 
                  (setf (aref v i) 
                        (float (funcall function i n) 0.0d0))))))))

(defun clip-in-window (x start end) (max start (min x end)))

(defun extract-window-into (vector start length destination)
  "Copy an extent of VECTOR to DESTINATION. Outside of its legal array
indices, VECTOR is considered to be zero."
  (assert (<= length (length destination)))
  (let ((start* (clip-in-window start 0 (length vector)))
        (end*   (clip-in-window (+ start length) 0 (length vector))))
    (unless (= length (- end* start*))
      (fill destination (coerce 0 (array-element-type destination))))
    (when (< -1 (- start* start) (length destination))
      (replace destination vector
               :start1 (- start* start)
               :end1 (+ (- start* start) (- end* start*))
               :start2 start*
               :end2 end*)))
  destination)

(defun extract-window 
    (vector start length &optional (element-type (array-element-type vector)))
  (extract-window-into 
   vector start length
   (make-array length
               :initial-element (coerce 0 element-type)
               :element-type element-type
               :adjustable nil
               :fill-pointer nil)))

(defun extract-centered-window-into (vector center size destination)
  "Extract a subsequence of SIZE from VECTOR, centered on OFFSET and
padding with zeros beyond the boundaries of the vector, storing it to
DESTINATION."
  (extract-window-into vector (- center (floor size 2)) size destination))

(defun extract-centered-window 
    (vector center size &optional (element-type (array-element-type vector)))
  "Extract a subsequence of SIZE from VECTOR, centered on CENTER and
padding with zeros beyond the edges of the vector."
  (extract-centered-window-into 
   vector center size
   (make-array size
               :initial-element (coerce 0 element-type)
               :element-type element-type
               :adjustable nil
               :fill-pointer nil)))

(defun convert-to-complex-sample-array (array)
  (let ((output (make-array (length array) 
                            :element-type 'complex-sample 
                            :adjustable nil
                            :fill-pointer nil)))
    (typecase array
      ((simple-array single-float 1)
       (loop for index from 0 below (length array)
             for x across array
             do (setf (aref output index) (complex (float x 0.0d0) 0.0d0))))
      ((simple-array double-float 1)
       (loop for index from 0 below (length array)
             do (setf (aref output index) (complex (aref array index) 0.0d0))))
      (t
       (loop for index from 0 below (length array)
             for x across array
             do (setf (aref output index) (complex (float (realpart x) 0.0d0) 0.0d0)))))
    output))

(defun windowed-fft (signal-vector center length &optional (window-fn 'hann))
  "Perform an FFT on the window of a signal, centered on the given
index, multiplied by a window generated by the chosen window function"
  (declare (type fixnum length)
           (optimize (speed 3)))
  (unless (power-of-two length)
    (error "FFT size ~D is not a power of two" length))
  (unless (typep signal-vector 'complex-sample-array)
    (setf signal-vector (convert-to-complex-sample-array signal-vector)))  
  (let* ((input-window (extract-centered-window signal-vector center length))
         (window (the (simple-array double-float 1)
                   (window-vector window-fn length))))
    (declare (type complex-sample-array input-window))
    (loop for index from 0 below length do
          (setf (aref input-window index)
                (* (aref input-window index)
                   (aref window index))))
    (fft input-window)))
