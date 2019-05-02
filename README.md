# Some notes

- There were some temporary, subroutine-local fluxes that I need to
  change to being in WaterFluxType, because now we need a separate
  version for each tracer. (I don't love this, but I can't see a clean
  way around it.)

# Some questions: organizing into subroutines

- How much should be split out into subroutines vs. staying in the
  top-level subroutine? Everything? Just state updates? Just tracer
  updates? Just state updates + tracer updates?
  
- For subroutines that operate on all tracers (or bulk + all tracers):
  should the loop over tracers be inside or outside the subroutine? Note
  that I have done something different for the allBrokenOut
  vs. tracersAndStatesBrokenOut options. I did what made the most sense
  to me for each of these, but it doesn't need to be the way it is for
  either of those. Also note: having the loop outside the subroutine
  makes it simpler to have associates for individual variables inside
  the subroutine. It's possible to have associates for individual
  variables otherwise, but this should probably be done via nested
  associates, which gets a bit more complex.
  
- Should we pass derived types into the subroutines or individual
  arrays? The former leads to simpler, more stable interfaces, but the
  latter does a better job at showing data flow.

- What (if anything) of the split-out things should be lumped together
  into a single subroutine, vs. each logically coherent thing being in
  its own subroutine? (For example: we could combine the update of
  snocan based on snow unloading with the computation of summed fluxes
  onto ground, because both operate over bulk + all tracers. But
  *should* we?)
  
# Reasons I was at least initially inclined to having things broken out

- I like the rule that routines that have interesting science shouldn't
  be complicated / obscured by tracer-related code. This is violated by
  the inline option, and to a lesser extent by the
  tracersAndStatesBrokenOut option.

- I like having routines either be coordination routines or calculation
  / sciencey routines, where the latter only operate on either bulk or a
  single tracer. This keeps these calculation routines easier to
  understand, and less prone to errors involving doing a calculation
  that is mixing different tracers or mixing the bulk with a tracer (as
  we are prone to for the ciso code). When you mix high-level
  coordination with lower-level science equations, I find it harder to
  understand: if you're interested in the high-level logic, you get
  distracted by the low-level equations, and vice versa

- Putting state updates in their own routine is a step closer to being
  able to separate physics from numerics (at least, I *think* that will
  help long-term)

# Question about calculating the tracer version of fluxes

For canopy interception and throughfall: Note that we compute a
fraction, then get the bulk flux by multiplying the source by this
fraction. We could also recast some of the other fluxes this way - e.g.,
I think this would be pretty straightforward for `qflx_snow_unload`, and
in principle, I think this would be possible for many fluxes.

The question is: for these fluxes, should we 

1. Keep the code as it's written now, where we calculate a bulk flux and
   then calculate a tracer flux based on the ratio of tracer source to
   bulk source

2. Recast the code to look like:

   - Calculate and save the fraction based on bulk quantities
   
   - In a loop over bulk and all tracers, do: `foo = frac * bar`
   
Overall, I lean towards (1), but I see arguments both ways.

Arguments for (1):

- It keeps the code simpler for bulk: fluxes are completely computed in
  a bulk-only section, rather than having the fraction computed in a
  bulk-only section then having the fluxes computed in a bulk-and-tracer
  section

- It keeps the method for setting tracer fluxes more consistent, because
  presumably we would *not* do method (2) for every flux in the model

Arguments for (2):

- It is less prone to the error of using the wrong source state in the
  tracer updates

- It may perform better in some cases
