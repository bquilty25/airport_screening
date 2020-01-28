# generate the current date and sink it into an MD file

sprintf("Last built at %s by B. Quilty, S. Clifford, S. Flasche, R. Eggo and other members of CMMID at LSHTM", 
        format(Sys.time(),
               "%d %b %Y at %H:%M:%S")) %>%
  write(., file = "date_stamp.md")
