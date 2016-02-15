
(defparameter *n* 500)
(defparameter *ci* .85)
(defparameter *bmax* .22)
(defparameter *bmin* .04)
(defparameter *tau* .13)

(defun transmission (u)
  (+ *bmax* (* u (- *bmin* *bmax*))))
(defun delta-infectives (i u)
  (let* ((s (- *n* i))
         (inew (round (/ (* (transmission u) i s) (1- *n*)))) 
         (snew (round (* *tau* i)))
         (splus (- (+ s snew) inew))  
         (iplus (- *n* splus)))
    (- iplus i)))

(defun delta-cost (i u &optional (sc 0))
  (let* ((s (- *n* i))
         (cu (/ *n* (+ 1 i)))
         (inew (round (/ (* (transmission u) i s) (1- *n*))))
         (c (+ (* *ci* inew) (* cu u))))
    (+ sc c)))

(defun leaf-cost (program &optional (infectives 100))
  (let ((results `((,infectives 0)))
        (cost 0))
    (dolist (u program (nreverse results))
      (incf cost (delta-cost infectives u))
      (incf infectives (delta-infectives infectives u))
      (push `(,infectives ,cost) results))))

(and 
  (< (abs (- (delta-cost 100 0) 15.3     )) 0.0001)
  (< (abs (- (delta-cost 100 1) 7.5004954)) 0.0001)
  (= (delta-infectives 100 0) 5)
  (= (delta-infectives 100 1) -10)
(= 2755 (apply #'+ (mapcar #'(lambda (tup) (car tup)) (leaf-cost '(1 1 1 1 1 1 1 1 1 1 0 0 0 1 1 0 0 0 1
             1 0 0 0 1 1 0 0 0 1 1 0 0 0 1 1 0 0 0
             1 1 0 0 0 1 1 0 0 0 1 1 0 0 0 1 1 0 0
             0 0 0 0)
           100))))
(<  (abs (- 17440.814  (apply #'+ (mapcar #'(lambda (tup) (cadr tup)) (leaf-cost '(1 1 1 1 1 1 1 1 1 1 0 0 0 1 1 0 0 0 1
             1 0 0 0 1 1 0 0 0 1 1 0 0 0 1 1 0 0 0
             1 1 0 0 0 1 1 0 0 0 1 1 0 0 0 1 1 0 0
             0 0 0 0)
           100))))) 0.0001))

(let ((best-cost 999999999)
      (best-infectives 0)
      (best-pgm nil)
      (leaves 0)
      (nodes 0))
  (defun update-best (cost infectives)
    (when (< cost best-cost)
      (setq best-cost cost)
      (setq best-infectives infectives)))
  (defun reset-best ()
    (setq best-cost 999999)
    (setq best-infectives 0)
    (setq best-pgm nil)
    (setq leaves 0)
    (setq nodes 0))
  (defun get-best-cost ()
    best-cost)
  (defun get-best-pgm ()
    best-pgm)
  (defun inc-leaves ()
    (incf leaves))
  (defun get-leaves ()
    leaves)
  (defun inc-nodes ()
    (incf nodes))
  (defun get-nodes ()
    nodes))

(defun get-cost (v)
  (first v))
(defun get-infectives (v)
  (second v))
(defun get-pgm (v)
  (third v))
(defun make-solution (cost infectives pgm)
  ;;(format t "~%make-solution ~a ~a ~a" cost infectives pgm)
  (list cost infectives pgm))

(defun dig (horizon infectives cost budget branches H)
  (let ((optimum (list 999999 0 nil))
        (current 0))
    (dotimes (i branches optimum)
      (let ((u (/ i (1- branches))))
        (when (< u budget)
          (let* ((delta-i (delta-infectives infectives u))
                (delta-c (delta-cost infectives u))
                (child (best-control
                        (1- horizon)
                        (+ infectives delta-i)
                        0
                        (- budget u)
                        branches
                        H)))
            (when (< (setq current (+ cost delta-c (get-cost child))) (get-cost optimum))
              (setq optimum (make-solution
                             current
                             (get-infectives child)
                             (cons u (get-pgm child)))))))))))
                
(defun best-control (horizon infectives cost budget branches &optional H)
  "Returns the final cost, infectives  and the corresponding control program."
  (if (zerop horizon)
      (progn
        (inc-leaves)
        (update-best cost infectives)
        (list cost infectives nil))
      (progn 
        ;; Not a leaf
        (inc-nodes)
        (when (null H)
          (setq H (make-hash-table :test 'equal)))
        (let* ((key `(,horizon ,infectives ,budget))
               (result (gethash key H)))
          (when (not result)
            (setq result (dig horizon infectives cost budget branches H))
            (setf (gethash key H) result))
          result))))

(defun main (&optional (horizon 60) (infectives 100) (budget 30)  (branches 2))
  (reset-best)
  (format t "~%Starting search~%Horizon: ~a~%Infectives: ~a~%Budget:~a" horizon infectives budget) 
  (let ((solution (best-control horizon infectives 0 budget branches)))
    (format t "~%Best cost: ~a ~%Best infectives: ~a~%Best pgm: ~a ~%Leaves: ~a~%Nodes:~A"
            (get-cost solution) (get-infectives solution) (get-pgm solution) (get-leaves) (get-nodes))))

(main)
