;; http://www.cs.cmu.edu/~mws/rpos.html
;; cos(arcsin(z))=sqrt(1-z^2)

(eval-when (:compile-toplevel)
 (setf asdf:*central-registry* (list "/home/martin/0505/mma/"))
 (require :gui))

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
	     (sin-theta (sqrt (- 1 (* z z))))
	     (sin-phi (sin phi))
	     (cos-phi (sqrt (- 1 (* sin-phi sin-phi)))))
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
				2 #+nil 3)
			 :element-type 'single-float))
	 (k 0))
    (dotimes (i na)
      (let ((x (aref a i 0))
	    (y (aref a i 1))
	   #+nil (z (aref a i 2)))
       (dotimes (j nb)
	 (setf (aref res k 0) (+ x (aref b j 0))
	       (aref res k 1) (+ y (aref b j 1))
	       ;(aref res (* i j) 2) (+ z (aref b j 2))
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

(loop for i below 100 collect
     (floor (* 20 (sqrt (+ (expt (aref *c* i 0) 2)
		    (expt (aref *c* i 1) 2))))))

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
		   sv #+nil(/ sv
		      r)
		   (/ (* 9 (sqrt sv)) (* r r)))))))))

#+nil
(time
 (let ((i 0))
   (loop for ratio in '(.98 .95 .8 .5 .01) do
	(let* ((theta-max #.(coerce (* (/ pi 180) 8) 'single-float))
	       (a (random-sphere-ring :theta-max theta-max
				      :theta-min 0s0 #+nil (* ratio theta-max) :n 2000))
	       (dz (abs (- (cos (* ratio theta-max)) ;; height of the ring cutout from sphere
			   (cos theta-max))))
	       (rmax (* 2.1 (sin theta-max)))
	       (w 130))
	  (output (project-and-bin (point-cloud-convolve a a)
				   :w w :rmax rmax)
		  ; :scale (/ (sphere-ring-area 1 dz))
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


(defparameter *b* (random-unit-disk :n 3000))
(defparameter *c* (point-cloud-convolve *b* *b*))
(defparameter *rot* 0)
(let ((a 'n))
  (defun draw-screen ()
    (gl:clear :color-buffer-bit)
    (gl:load-identity)
    (gl:point-size 3)
    (gl:blend-func :src-alpha :one-minus-src-alpha)
    
    (gl:scale 60 60 60)
    (incf *rot*)
    (gl:translate 2 2 0)
    (gl:rotate *rot* 0 0 1)
   
    (gl:enable :blend)
    (gl:color 1 1 1 .3)
    (gl:with-primitive :points
      
      (dotimes (i (array-dimension *b* 0))
	(gl:vertex (aref *b* i 0)
		   (aref *b* i 1)
		   #+nil(aref *c* i 2))))))

#+nil
(let ((x 700)
      (y 100))
  (sb-thread:make-thread
   #'(lambda ()
       (gui:with-gui (300 300 x y)
         (draw-screen)))
   :name "display"))
