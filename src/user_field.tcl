
if { [ info exists user_data ] == 0 } { 

  puts "#  "
  puts "# user_field.tcl - from pdblovofit "
  puts "#  "
  puts "# ERROR:  "
  puts "#  "
  puts "# The variable \"user_data\" must be defined first, containing"
  puts "# the name of the pdb file created by pdblovofit."
  puts "#  "
  puts "# Use, for example:  set user_data ./align.pdb "
  puts "#  "

} else {

  set numframes [ molinfo top get numframes ]
  set numatoms [ molinfo top get numatoms ]
  
  set user_data [ open $user_data r ] 
  set user_data [ read $user_data ]
  set user_data [ split $user_data "\n" ]
  set iframe 0
  set iall 0
  foreach line $user_data {
  
    if { [ string range $line 0 3 ] == "ATOM" |
         [ string range $line 0 5 ] == "HETATM" } {       
  
      if { $iall == 0 } {
        animate goto $iframe
        puts "Setting User fields for frame: $iframe"
      }
      set user  [ string range $line 54 59 ]
      set user2 [ string range $line 60 65 ]
  
      set atomselect [ atomselect top "index $iall" ]
      $atomselect set user $user 
      $atomselect set user2 $user2
      $atomselect delete
  
      incr iall
      if { $iall == $numatoms } { 
        incr iframe
        set iall 0
      } 
    } 
  }

}

