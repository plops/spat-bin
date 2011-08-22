;; http://www.cs.cmu.edu/~mws/rpos.html
;; cos(arcsin(z))=sqrt(1-z^2)

(declaim (optimize (speed 3) (safety 1) (debug 1)))

(defun random-sphere-ring (&key (theta-max #.(coerce (* 67 (/ 2 pi 180))
						     'single-float))
			   (theta-min 0s0) (n 100))
  (declare (type single-float theta-min theta-max)
	   (type fixnum n))
  (let* ((res (make-array (list n 3) :element-type 'single-float))
	 (zmin (cos (max theta-max
			 theta-min)))
	 (zmax (cos (min theta-max
			 theta-min)))
	 (dz (- zmax zmin)))
    (assert (<= zmax 1s0))
    (dotimes (i n)
      (let* ((z (+ zmin (random dz)))
	     (phi (random (coerce (* 2 pi) 'single-float)))
	     (cos-theta (sqrt (- 1 (* z z))))
	     (sin-phi (sin phi))
	     (cos-phi (sqrt (- 1 (* sin-phi sin-phi)))))
	(declare (type (single-float 0s0 1s0) z))
	(setf (aref res i 0) (* cos-theta cos-phi)
	      (aref res i 1) (* cos-theta sin-phi)
	      (aref res i 2) z)))
    res))


#+nil
(random-sphere-ring)

(deftype point-cloud ()
  `(simple-array single-float 2))

(defun point-cloud-convolve (a b)
  (declare (type point-cloud a b)
	   (values point-cloud &optional))
  (let* ((na (array-dimension a 0))
	 (nb (array-dimension b 0))
	 (res (make-array (list (* na nb) 2 #+nil 3)
			 :element-type 'single-float)))
    (dotimes (i na)
      (let ((x (aref a i 0))
	    (y (aref a i 1))
	   #+nil (z (aref a i 2)))
       (dotimes (j nb)
	 (setf (aref res (* i j) 0) (+ x (aref b j 0))
	       (aref res (* i j) 1) (+ y (aref b j 1))
	       ;(aref res (* i j) 2) (+ z (aref b j 2))
	       ))))
    res))

#+nil
(let ((a (random-sphere-ring)))
  (point-cloud-convolve a a))

;; http://en.wikipedia.org/wiki/Spherical_cap
(defun sphere-ring-area (r dz)
  (* #.(coerce (* 2 pi) 'single-float) r dz))

(defun project-and-bin (a &key (w 128) (rmax 2s0))
  (let ((hist (make-array w :element-type '(unsigned-byte 64))))
    (declare (type (simple-array (unsigned-byte 64) 1) hist))
    (dotimes (i (array-dimension a 0))
      (incf (aref hist
		  (min (1- w) (round (* (/ w rmax)
					(sqrt
					 (+ (expt (aref a i 0) 2)
					    (expt (aref a i 1) 2)))))))))
    
    (loop for i below w collect
	 (list (* i (/ rmax w)) (aref hist i)))))

(defun output (a &key (scale 1s0) (append nil))
  (with-open-file (s "/dev/shm/rad.dat" :direction :output
		     :if-exists (if append
				    :append
				    :supersede)
		     :if-does-not-exist :create) 
    (format s "~%")
    (dotimes (i (length a))
      (unless (= i 0)
	(destructuring-bind (r v) (elt a i)
	  (let ((sv (* scale v)))
	   (format s "~f ~f ~f~%"
		   r
		   (/ sv
		      (* r r))
		   (/ (* 9 (sqrt sv)) (* r r)))))))))

#+nil
(time
 (let ((i 0))
   (loop for ratio in '(.98 .95 .8 .5 .01) do
	(let* ((theta-max #.(coerce (* (/ pi 180) 68) 'single-float))
	       (a (random-sphere-ring :theta-max theta-max
				      :theta-min (* ratio theta-max) :n 2000))
	       (dz (abs (- (cos (* ratio theta-max)) ;; height of the ring cutout from sphere
			   (cos theta-max))))
	       (rmax (* 2.1 (sin theta-max)))
	       (w 130))
	  (output (project-and-bin (point-cloud-convolve a a)
				   :w w :rmax rmax)
		  :scale (/ (sphere-ring-area 1 dz))
		  :append (/= i 0)))
	(incf i))))
