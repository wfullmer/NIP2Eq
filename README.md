# NIP2Eq


# old note
This was an old note that I left to myself when I wrote this code back in 2012 or so:

need to add: if(bc = periodic) then call bndry     affter the source call
  then wrap all the "update" terms into a subroutine
  need to apply bndry after every exact call (see exact subroutine, all calls need to be removed)
  put all subroutines into a seperate file that will get cat'd into a large otherwise empty module
