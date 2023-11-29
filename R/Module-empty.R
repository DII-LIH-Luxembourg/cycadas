
Blah_UI <- function(blah) {
  ns <- NS(blah)

}

Blah_Server <- function(blah) {
  moduleServer(
    id = blah,
    module = function(input, output, session) {
 
    })}