;; http://www.cs.cmu.edu/~mws/rpos.html
;; cos(arcsin(z))=sqrt(1-z^2)

(defun random-sphere-ring (&key (theta-max #.(coerce (* 67 (/ 2 pi 180))
						     'single-float))
			   (theta-min 0) (n 100))
  (let* ((res (make-array (list n 3) :element-type 'single-float))
	 (zmin (cos (max theta-max
			 theta-min)))
	 (zmax (cos (min theta-max
			 theta-min)))
	 (dz (- zmax zmin)))
    (dotimes (i n)
      (let* ((z (+ zmin (random dz)))
	     (phi (random (coerce (* 2 pi) 'single-float)))
	     (cos-theta (sqrt (- 1 (* z z))))
	     (sin-phi (sin phi))
	     (cos-phi (sqrt (- 1 (* sin-phi sin-phi)))))
	(setf (aref res i 0) (* cos-theta cos-phi)
	      (aref res i 1) (* cos-theta sin-phi)
	      (aref res i 2) z)))
    res))


#+nil
(random-sphere-ring)

(defun point-cloud-convolve (a b)
  (let* ((na (array-dimension a 0))
	 (nb (array-dimension b 0))
	 (res (make-array (list (* na nb) 3)
			 :element-type 'single-float)))
    (dotimes (i na)
      (let ((x (aref a i 0))
	    (y (aref a i 1))
	    (z (aref a i 2)))
       (dotimes (j nb)
	 (setf (aref res (* i j) 0) (+ x (aref b j 0))
	       (aref res (* i j) 1) (+ y (aref b j 1))
	       (aref res (* i j) 2) (+ z (aref b j 2))))))
    res))

#+nil
(let ((a (random-sphere-ring)))
  (point-cloud-convolve a a))

(defun project-and-bin (a &key
			(w 128) (h 128)
			(ox (/ w 2s0))
			(oy (/ h 2s0))
			(sx (/ w 2s0))
			(sy (/ h 2s0)))
  (let ((res (make-array (list h w)
			 :element-type '(unsigned-byte 64))))
    (dotimes (i (array-dimension a 0))
      (incf (aref res 
		  (min (1- w) (max 0 (round (+ ox (* sx (aref a i 0))))))
		  (min (1- h) (max 0 (round (+ oy (* sy (aref a i 1)))))))))
    (let* ((ma (reduce #'max (sb-ext:array-storage-vector res)))
	   (s (/ 255s0 ma))
	   (img (make-array (list h w) :element-type 'unsigned-byte)))
      (dotimes (j h)
	(dotimes (i w)
	  (setf (aref img j i) (floor (* s (aref res j i))))))
      img)))

#+nil
(let ((a (random-sphere-ring :theta-min .1 :n 1000)))
  (project-and-bin
   (point-cloud-convolve a a)
   :w 16 :h 16))