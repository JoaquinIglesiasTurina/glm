(defproject stats-projects "0.1.0-SNAPSHOT"
  :description "FIXME: write description"
  :url "http://example.com/FIXME"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[org.clojure/clojure "1.8.0"]
                 [incanter/incanter-core "1.9.1"]
                 [incanter/incanter-core "1.9.1"]                 
                 [incanter/incanter-core "1.9.1"]
                 [incanter/incanter-charts "1.9.1"]
                 [incanter/incanter-excel "1.9.1"]]
  :main ^:skip-aot stats-projects.core
  :target-path "target/%s"
  :profiles {:uberjar {:aot :all}})

