(ns stats-projects.core
  (:require [clojure.java.io :as io]
            [incanter.core :as i]
            [incanter.stats :as istats]                        
            [incanter.optimize :as iopt]
            [incanter.excel :as xls]))
;; TODO: test aic, test deviance analysis
;; utility functions
(defn apply-to-all
  "applies a function to all elements of a nx1 matrix (math vector)"
  [f vector]
  (reduce (fn [a-seq elem]
            (let [to-apply (first elem)
                  applied (f to-apply)]
              (conj a-seq applied)))
          [] vector))

(defn discretize
  "turns real value to 1 or 0. Used for testing"
  [number]
  (if (> number 0.5)
    1
    0))

(defn to-matrix
  "turns vector into nx1 matrix. makes it easier to use the optimizer."
  [a-seq]
  (if (i/matrix? a-seq)
    a-seq
    (i/matrix a-seq 1)))

;; define families' link and log-likelihood functions
; binomial family
(defn logit
  "Logit function. Takes fitted data"
  [k]
  (let [expk (i/exp k)]
    (i/div expk (i/plus 1 expk))))

(defn probit
  "Probit function, takes fitted data"
  [k]
  (i/matrix (apply-to-all istats/cdf-normal k) 1))

(defn log-binomial
  "log-binomial density function takes link, 
  data (dependent and independent) and parameters"
  [link y x theta]
  (defn helper [y p]
    (cond (== y 1) (i/log p)
          (== y 0) (i/log (- 1 p))
          :else 0))
  (let [theta' (to-matrix theta)
        y-star-fitted (i/mmult x theta')
        p' (apply-to-all identity (link y-star-fitted))
        y' (apply-to-all identity y)]
    (/ (i/sum (map helper y' p')) (i/nrow x))))

;; multimethods for link gradient and variance function
; gradients
(defmulti gradient :link)

(defmethod gradient logit
  [parametrization]
  (let [theta' (to-matrix (:beta parametrization))
        x (:x parametrization)
        y-star (i/mmult x theta')
        p (logit y-star)]
    (i/mult -1 (i/bind-columns p p) x)))

(defmethod gradient probit
  [parametrization]
  (let [theta' (to-matrix (:beta parametrization))
        x (:x parametrization)
        y-star (i/mmult x theta')]
    (i/matrix (apply-to-all istats/pdf-normal y-star) 1)))

;variance
(defmulti cov-matrix :family)

(defmethod cov-matrix log-binomial
  [parametrization]
  (let [link (:link parametrization)
        theta' (to-matrix (:beta parametrization))
        x (:x parametrization)
        p (link (i/mmult x theta'))
        p' (gradient parametrization)
        w (i/mult p (i/minus 1 p))]
    (i/solve (i/mmult (i/trans p') (i/mult (i/bind-columns w w) p')))))

;; families map
(def families {:binomial
               {:function log-binomial
                :links {:logit logit :probit probit}}})
;; generalized linear model
(defn glm [y x start & {:keys [family link]
                        :or {family log-binomial link logit}}]
  ;log-likelihood helper function, reduces verbosity
  (defn log-lik [x theta]
    (family link y x theta))
  ;optimization helper function, reduces verbosity on deviance analysis
  (defn max-lik [x start]
    (let [objective (fn [theta] (log-lik x theta))]
      (:value (iopt/maximize objective start))))
  (let [n (i/nrow x) ;number of obvservations
        ;xs for deviance analysis
        [ones id] [(i/matrix (repeat n 1) 1) (i/identity-matrix n)]        
        beta-hat (max-lik x start)
        log-model (log-lik x beta-hat)
        ;betas and log-like for deviance analysis
        [beta-null beta-full] [(max-lik ones [0]) (max-lik id (repeat n 0))]
        [log-null log-full] [(log-lik ones beta-null) (log-lik id beta-full)]
        [null-deviance residual-deviance]
        (mapv (fn [log] (i/mult 2 (- log-full log))) [log-null log-model])
        ;computes variance
        parametrization {:x x :y y :beta beta-hat :link link
                         :family family} ;parametrization map for variance
        variance (i/diag (cov-matrix parametrization)) ;compute variance
        ;akaike information criterion
        aic (* 2 (- n (family link y x beta-hat)))
        ]
    {:beta beta-hat :variance variance :aic aic
     :null-deviance null-deviance :residual-deviance residual-deviance
     }))

;; simulations and tests
; covariates
(def x (istats/sample-mvn 10 :sigma (i/diag [1 2])))
; parameters
(def beta (i/matrix [0.5 -0.3] 1))
; latent variable
(def y-star (i/plus (i/mmult x beta)
                    (istats/sample-mvn (i/nrow x) :sigma (i/diag  [1]))))
; obvserved variable
(def y (i/matrix (apply-to-all discretize (logit y-star)) 1))
; parametrization for variance multi
(def parametrization {:x x :y y :beta beta :link logit :family log-binomial})
; test optimization
(def optzed (glm y x [0 0]))
(println optzed)


