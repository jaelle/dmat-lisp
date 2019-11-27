; The following is adapted from the DMA Toolbox (https://ppw.kuleuven.be/okp/software/dmat/)

(ql:quickload "cl-randist")

(defvar *outfile* "simuldiff_results.csv")

(defun start (par n)
 
  (with-open-file (out *outfile* :direction :output :if-exists :supersede :if-does-not-exist :create)
    (format out "run,RT,correct,~%"))
    
  (loop
     for i from 1 to n
     do
       (simuldif par i)))
  


(defun simuldif (par i)
	; must have 7 parameters
  (if (not (equal (length par) 7))
      nil)

					; drift rate must be between [-0.5,0.5]
  (let ((nu (nth 6 par)))
    (if (or
	 (< nu -0.5)
	 (> nu 0.5))
	nil))
       
  (let* ((a (nth 0 par))
	 (z (nth 3 par))
	 (sz (nth 4 par))
	 
	 (zz (+ (- z (/ sz 2) (* sz (randist:random-uniform)))))
	 (Aupper (- a zz))
	 (Alower (- 0 zz))
	 (radius (min (abs Aupper) (abs Alower)))
	 (results (wiener-process par 0 Aupper Alower radius 0)))

	(print results)

    (with-open-file (out *outfile* :direction :output :if-exists :append :if-does-not-exist :create)
      (format out "~A,~,6f,~A,~%"
	      i
	      (nth 0 results)
	      (nth 1 results)))))

(defun wiener-process(par startpos Aupper Alower radius totaltime)

  (let* ((a (nth 0 par))
	 (Ter (nth 1 par))
	 (eta (nth 2 par))
	 (z (nth 3 par))
	 (sz (nth 4 par))
	 (st (nth 5 par))
	 (nu (nth 6 par))
	 (tau 0.1)
	 (sigma (/ (* tau tau) 2))
	 
	 (delta 0.0000001)
	 
	 (r1 (randist:random-normal))
	 ;(r1 0.5)
	 (mu (+ nu (* r1 eta)))
	 
	 (dlambda (+ (/ (* 0.25 (expt mu 2)) sigma) (/ (* 0.25 (* sigma (expt pi 2))) (expt radius 2))))

	 (F1 (/ (* sigma pi) (* radius mu)))
	 (F (/ (expt F1 2) (+ 1 (expt F1 2))))
	 
	 (prob1 (exp (/ (* radius mu) sigma)))
	 (prob (/ prob1 (+ 1 prob1)))
	 (randval (randist:random-uniform))
	 ;(randval 0.5)

	 (dir1 (if (< randval prob)
		   1
		   -1))
	 (dir2 (+ (* dir1 radius) startpos))
	 
	 (totaltime2 (+ totaltime (simultime F dlambda)))
	 ;(ndrt (+ (- Ter (/ st 2) (* st 0.5))))
	 (ndrt (+ (- Ter (/ st 2) (* st (randist:random-uniform)))))
	 (finaltime (+ totaltime2 ndrt)))
	 
	 (cond ((> (+ dir2 delta) Aupper)
		(list finaltime 1))
	       ((< (- dir2 delta) Alower)
		(list finaltime 0))
	       (t
		(wiener-process par dir2 Aupper Alower (min (abs (- Aupper dir2)) (abs (- Alower dir2))) totaltime2)))))

(defun simultime (F dlambda)

  (let ((s1 (simultime1 0 0 -1 F)))

    (/ (abs (log s1)) dlambda)))


(defun simultime1 (s1 s2 l F)
  (cond ((> s2 l) ; recursive case

      (let* ((new_s2 (randist:random-uniform))
	     (new_s1 (randist:random-uniform))
      ;(let* ((new_s2 0.5)
	     ;(new_s1 0.5)
	     (tnew (simultime2 0 0 0 new_s1 F))
	     (new_l (+ 1 (* (expt new_s1 (- 0 F)) tnew))))

	  (simultime1 new_s1 new_s2 new_l F)))

	  (t s1)))

(defun simultime2 (tnew told uu s1 F)
  (cond ((or
       (> (abs (- tnew told)) 0.0000001)
       (equal uu 0)) ;recursive case

      (let* ((told2 tnew)
	     (new_uu (+ uu 1))
	     (tnew2_pt1 (+ (* 2 new_uu) 1))
	     (tnew2_pt2 (expt -1 new_uu))
	     (tnew2_pt3 (expt s1 (* F (expt (+ (* 2 new_uu) 1) 2))))
	     
	     (tnew2 (+ told2 (* (* tnew2_pt1 tnew2_pt2) tnew2_pt3))))
	;s(break)
	(simultime2 tnew2 told2 new_uu s1 F)))
	
	(t tnew)))
	    
      

(defun set-eta (eta)
  
					; drift rate standard deviation must be <= 0.3  
       (if (> eta 0.3)
	   nil)

       
       (if (equal eta 0)
	   1e-16))

