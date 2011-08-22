;; http://www.cs.cmu.edu/~mws/rpos.html
;; cos(arcsin(z))=sqrt(1-z^2)

(eval-when (:compile-toplevel)
 (setf asdf:*central-registry* (list "/home/martin/0505/mma/"))
 (require :gui))

#+nil(declaim (optimize (speed 3) (safety 1) (debug 1)))

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
	     (phi (random #.(coerce (* 2 pi) 'single-float)))
	     (sin-theta (sin (acos z)))
	     (sin-phi (sin phi))
	     (cos-phi (cos phi)))
	(declare (type (single-float 0s0 1s0) z))
	(setf (aref res i 0) (* sin-theta cos-phi)
	      (aref res i 1) (* sin-theta sin-phi)
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
	 (res (make-array (list (* na nb) 
				3)
			 :element-type 'single-float))
	 (k 0))
    (dotimes (i na)
      (let ((x (aref a i 0))
	    (y (aref a i 1))
	    (z (aref a i 2)))
       (dotimes (j nb)
	 (setf (aref res k 0) (+ x (aref b j 0))
	       (aref res k 1) (+ y (aref b j 1))
	       (aref res k 2) (+ z (aref b j 2))
	       )
	 (incf k))))
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
      (let ((r (sqrt
		(+ (expt (aref a i 0) 2)
		   (expt (aref a i 1) 2)))))
       (incf (aref hist
		   (max 0 (min (1- w) (floor (* (/ w rmax)
						r))))))))
    
    (loop for i below w collect
	 (list (* i (/ rmax w)) 
	       (aref hist i)))))



#+nil
(mapcar #'(lambda (x) (second x))
	(project-and-bin *b*
			 :rmax 1.5s0
			 :w 15))


(defun project-and-bin-2 (a &key (w 128) (h 128) (xmax 4s0) (ymax 4s0))
  (let ((hist (make-array (list h w) :element-type '(unsigned-byte 64))))
    (declare (type (simple-array (unsigned-byte 64) 2) hist))
    (dotimes (i (array-dimension a 0))
      (incf (aref hist 
		  (max 0 (min (1- h) (round (* (/ h ymax) (+ 3 (aref a i 1))))))
		  (max 0 (min (1- w) (round (* (/ w xmax) (+ 3 (aref a i 0)))))))))
    hist))

#+nil
(project-and-bin-2 *c* :w 16 :h 16 :xmax 6s0 :ymax 6s0)



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
	(let* ((theta-max #.(coerce (* (/ pi 180) 8) 'single-float))
	       (a (random-sphere-ring :theta-max theta-max
				      :theta-min (* ratio theta-max) :n 1000))
	       (dz (abs (- (cos (* ratio theta-max)) ;; height of the ring cutout from sphere
			   (cos theta-max))))
	       (rmax (* 2.1 (sin theta-max)))
	       (w 130))
	  (output (project-and-bin (point-cloud-convolve a a)
				   :w w :rmax rmax)
		  :scale (/ (sphere-ring-area 1 dz))
		  :append (/= i 0)))
	(incf i))))


(defun random-unit-disk (&key (n 100))
  (let ((i 0)
	(res (make-array (list n 3) :element-type 'single-float)))
    (loop while (< i n) do
	 (let* ((x (- (random 2s0) 1s0))
		(y (- (random 2s0) 1s0))
		(r2 (+ (* x x) (* y y))))
	   (when (< r2 1s0)
	     (setf (aref res i 0) x
		   (aref res i 1) y)
	     (incf i))))
    res))

#+nil
(time
 (let ((a (random-unit-disk :n 100)))
   (output
    (project-and-bin (point-cloud-convolve a a)
		     :w 128
		     :rmax 2.1s0))))

(defun invert-point-cloud (a)
  (let ((b (make-array (array-dimensions a)
		       :element-type 'single-float)))
    (dotimes (i (array-dimension a 0))
      (setf (aref b i 0) (- (aref a i 0))
	    (aref b i 1) (- (aref a i 1))
	    (aref b i 2) (- (aref a i 2))))
    b))
#+nil
(let ((theta-max (* #.(coerce (/ pi 180) 'single-float)
		    70)))
 (defparameter *b* (random-sphere-ring 
		    :theta-max theta-max
		    :theta-min (* .1 theta-max) :n 500)))
(defparameter *b* (random-unit-disk :n 2000))
(defparameter *c* (point-cloud-convolve *b* (invert-point-cloud *b*)))
(defparameter *rot* 0)
(let ((a 'n))
  (defun draw-screen ()
    (let ((w 500)
	  (h 500))
      (gl:load-identity)
      (gl:viewport 0 0 w h)
      (gl:matrix-mode :projection)
      (gl:load-identity)
      (progn (glu:perspective 40 (/ w h) 3 100)
	     (glu:look-at 20 30 -5
			  0 0 0
			  0 0 1))
      (gl:matrix-mode :modelview)
      (gl:load-identity))

    (gl:clear :color-buffer-bit :depth-buffer-bit)
    (gl:load-identity)
    (gl:point-size 3)
    (gl:blend-func :src-alpha :one-minus-src-alpha)
    
    (let ((s 5))
      (gl:scale s s s))
    (incf *rot*)
    
	

    (gl:rotate 10 1 0 0)
    (gl:rotate *rot* 0 0 1)

    
   
    (gl:enable :blend)
    (gl:disable :depth-test)
    
    (gl:with-primitive :lines
      (gl:color 1 0 0) (gl:vertex 0 0 0) (gl:vertex 1 0 0)
      (gl:color 0 1 0) (gl:vertex 0 0 0) (gl:vertex 0 1 0)
      (gl:color 0 0 1) (gl:vertex 0 0 0) (gl:vertex 0 0 1))
    (gl:color 1 1 1 .01)
    (gl:with-primitive :points
      (let ((e *c*))
       (dotimes (i (array-dimension e 0))
	 (gl:vertex (aref e i 0)
		    (aref e i 1)
		    (aref e i 2)))))))

#+nil
(let ((x 700)
      (y 100))
  (sb-thread:make-thread
   #'(lambda ()
       (gui:with-gui (500 500 x y)
         (draw-screen)))
   :name "display"))
