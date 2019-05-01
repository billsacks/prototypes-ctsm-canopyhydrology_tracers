# Some notes

- There were some temporary, subroutine-local fluxes that I need to
  change to being in WaterFluxType, because now we need a separate
  version for each tracer. (I don't love this, but I can't see a clean
  way around it.)

# Some questions

- How much should be split out into subroutines vs. staying in the
  top-level subroutine?
  
- What (if anything) of the split-out things should be lumped together
  into a single subroutine, vs. each logically coherent thing being in
  its own subroutine? (For example: we could combine the update of
  snocan based on snow unloading with the computation of summed fluxes
  onto ground, because both operate over bulk + all tracers. But
  *should* we?)
  
- For subroutines that operate on all tracers (or bulk + all tracers):
  should the loop over tracers be inside or outside the subroutine?
  
- Should we pass derived types into the subroutines or individual
  arrays? The former leads to simpler, more stable interfaces, but the
  latter does a better job at showing data flow.

