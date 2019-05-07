# Notes about this document

Most thoughts here are from Bill Sacks. Comments labeled "From
discussion" refer to a discussion (2019-05-02) between Bill and: Mike
Barlage, Erik Kluzek, Negin Sobhani, Sean Swenson and Mariana
Vertenstein.

# Some notes

- There were some temporary, subroutine-local fluxes that I need to
  change to being in WaterFluxType, because now we need a separate
  version for each tracer. (I don't love this, but I can't see a clean
  way around it.)
  
  **From discussion:** People are okay with this.
  
  Update (2019-05-07) See some more notes on this below: today's notes,
  under the heading, "Need to change some local, temporary variables to
  be part of a water type?"

# Some questions: organizing into subroutines

- How much should be split out into subroutines vs. staying in the
  top-level subroutine? Everything? Just state updates? Just tracer
  updates? Just state updates + tracer updates?
  
  **From discussion:** People are happy with everything broken
  out. (Note: we examined the inline and everything-broken-out
  options. We did not examine intermediate solutions because people were
  happy with the everything-broken-out option.) One reason people like
  this is because, if the names of subroutines are good, then you can
  easily get a high-level view of what the code is doing. Also, this
  lends itself better to unit tests.
  
- What (if anything) of the split-out things should be lumped together
  into a single subroutine, vs. each logically coherent thing being in
  its own subroutine? For example: we could combine the update of
  snocan based on snow unloading with the computation of summed fluxes
  onto ground, because both operate over bulk + all tracers. But
  *should* we? Similarly, we could combine the setting of tracer fluxes
  with the related state updates.
  
  **From discussion:** Let's start with what we have now, with more
  subroutines. We could always combine them more later.
  
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
  
  **From discussion:** No opinions. Update: based on the answer to the
  below question, I'll keep the loops outside the routines.
  
- Should we pass derived types into the subroutines or individual
  arrays? The former leads to simpler, more stable interfaces, but the
  latter does a better job at showing data flow. However, the latter is
  only possible if the loop over tracers is outside the call. We could
  also mix & match, passing individual arrays for stuff that works on
  bulk as well as tracers, but just passing `water_inst` for things that
  are tracer-only: this way you could at least see the data flow for
  bulk. In the solutions where we inline the bulk flux calculations:
  Passing individual arrays may have the most benefit for subroutines
  that involve bulk, but there is little enough work in those routines
  that, by the time we have the loop over instances and the
  many-argument calls, I think I'd prefer to just inline this
  code. However, passing individual arrays (at least for routines that
  operate on bulk water or bulk+tracers) could be a compromise for the
  allBrokenOut solution, allowing you to at least see data flow in that
  solution.
  
  **From discussion:** Pass individual arrays for anything that operates
  on bulk for (1) better ability to trace data flow, and (2) ability to
  declare variable intent. For things just calculating tracer fluxes, no
  strong feelings.

# Reasons I was at least initially inclined to having things broken out

- I like the rule that routines that have interesting science shouldn't
  be complicated / obscured by tracer-related code. This is violated by
  the inline option, to a lesser extent by the tracersBrokenOut option,
  and to an even lesser (but non-zero) extent by the
  tracersAndStatesBrokenOut option.

- I like having routines either be coordination routines or calculation
  / sciencey routines, where the latter only operate on either bulk or a
  single tracer. This keeps these calculation routines easier to
  understand, and less prone to errors involving doing a calculation
  that is mixing different tracers or mixing the bulk with a tracer (as
  we are prone to for the ciso code). When you mix high-level
  coordination with lower-level science equations, I find it harder to
  understand: if you're interested in the high-level logic, you get
  distracted by the low-level equations, and vice versa.

- Putting state updates in their own routine is a step closer to being
  able to separate physics from numerics (at least, I *think* that will
  help long-term).
  
- I worry that having everything inlined will lead to more degradation
  of the code over time: sections will grow and start to include things
  that they shouldn't, rather than keeping individual pieces more
  logically separated. (e.g., someone will put something in the state
  update section that doesn't belong there.)

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

**From discussion:** No real opinions

# Other notes from discussion

## Short-circuiting some state updates?

Sean: Can we short-circuit some state updates? In particular: For canopy
excess, can we avoid doing the update of liqcan, instead calculating the
excess and putting it in a runoff flux without ever updating liqcan to
be greater than the holding capacity?

Mike points out that this would have different implications for the
tracer concentration: the current scheme mixes together the inputs with
the existing state, whereas Sean's suggestion doesn't do this mixing
before the runoff.

We like Sean's idea for bulk, but aren't sure if that's the right thing
to do with respect to tracers. Let's come back to this....

## Naming

Rename PrecipInputs to make it clear that irrigation is added?

## Avoiding use statements

Let's pass col directly rather than using it.

## Comments about whether a subroutine is a state update vs. flux calculation

Comment whether something is a state update?

Actually, use prefixes on subroutine name to make this more clear.

# Notes after discussion

## (2019-05-03)

I think passing individual arguments greatly enhances the utility of the
top-level coordination routine. However, it does lead to a lot more
typing: (not counting the actual use in the science code) each variable
needs to be typed 5 times (2 in the call, 1 in the argument list, 1 in
the argument declaration, 1 in the shr_assert), rather than 2 times for
the alternative (2 times in the associate statement within a routine). I
favor read-time convenience over write-time convenience, so I think this
is a win, but what do others think?

I have used associate statements for 'b' and 'w'. The former lets us
access bulk and the latter lets us access the current bulk-or-tracer. I
normally prefer longer names, but here I went for conciseness to shorten
already-long lines. I'm happy to lengthen these if others want.

I lean slightly towards the solution that passes individual arrays for
all routines (not just those that operate on bulk). This adds long
argument lists for the tracer-only routines, but keeps the tracer-only
routines consistent with the bulk-only or bulk-and-tracer routines. In
addition to the consistency, I like this because I think this will make
it more likely for the sources (for state updates & tracer updates) to
stay in sync: I think it will be easier for someone to notice that they
need to change the argument list for the tracer flux update routine.

Note that, for the argument names to the tracer flux updating routine,
I'm using prefixes of `bulk_` and `trac_`. Using prefixes makes it more
obvious if you accidentally use a tracer variable where bulk should go,
or vice versa. Using `trac_` (rather than `tracer_`) makes the tracer
variables the same number of characters as the bulk, letting you see
that things are consistent at a glance.

## (2019-05-07)

### Need to change some local, temporary variables to be part of a water type?

#### Initial thoughts

There are some local, temporary variables
(`qflx_liq_above_canopy_patch`, `forc_snow_downscaled_patch`, and maybe
others) that are now part of a water type, and so fill up memory
throughout the run.

I can see how we could avoid this by computing these variables at the
last moment that we need them. This would be similar to what we do in
Wateratm2lndType.F90: SetOneDownscaledTracer (see also the initial
comment in https://github.com/ESCOMP/ctsm/issues/487): we set the
col-level variable equal to the gridcell-level variable at the last
possible moment. If we truly waited until the last possible moment,
though, we'd end up recomputing them repeatedly for the bulk. We could
store the bulk as a subroutine-local variable in the coordination
routine (just calculating the tracer versions of these two variables in
`TracerFlux_CanopyInterceptionAndThroughfall`); this would avoid
repeated calculation of the bulk quantity, but it would break symmetry
between the bulk and tracer in that tracer routine. In addition, it
would make it harder if we eventually want to defer these tracer
calculations until later in the driver loop, since the interface depends
on pre-calculated variables for the bulk which would no longer be
available. (Though we could get around that last point by doing a
compromise solution, where store the bulk quantity in the bulk-only
type, but calculate the tracer version locally.) Finally, I feel like
these various alternatives are generally more complicated and
confusing - e.g., if the bulk and tracer quantities are computed in
different places; this could also make it more likely that the tracer
and bulk calculations get out of sync.

So, at least for now, I'm going to stick with storing these variables in
the derived types. But I could imagine revisiting this, particularly if
we're running into this sort of thing a lot.

#### Updated thoughts

Upon further reflection, I see a way I like: First, just calculate these
quantities for the bulk. Then have the coordination routine calculate
these quantities for each tracer just before calling the
`TracerFlux_CanopyInterceptionAndThroughfall` routine (rather than
calling `SumFlux_TopOfCanopyInputs` from within
`TracerFlux_CanopyInterceptionAndThroughfall`, which would lead to the
asymmetry that I didn't like). This allows us to keep variables local
that should be local (which is good for code understandability as well
as memory use and performance), and keeps things fairly symmetrical
between the bulk and tracers.

This won't work if we ever move the tracer calculation outside of this
coordination routine, but (a) we can cross this bridge when we come to
it, and (b) maybe we'll want to keep the tracer calculations in this
coordination routine long-term anyway (e.g., maybe we'll want to
consolidate the current 3 TracerFlux routines into 1, but still keep
that 1 called from this coordination routine).

One downside of this is that the tracer quantities will always have
memory on the stack, even for runs without tracers. But I don't think
this is a huge deal. Another downside is that there is an additional
call to `SumFlux_TopOfCanopyInputs`, and it's slightly less obvious that
the tracers and bulk are doing the same thing.

Overall, I don't feel strongly about this, though, and could be
convinced to go back to the earlier idea of moving all of these
could-be-local variables into types.

However: I probably will NOT apply this to fluxes that appear in state
updates, because long-term, we may want to move these state updates to
some higher level.

What about `qflx_through_snow` and `qflx_through_rain`? These could
possibly remain local, since they aren't directly involved in state
updates: they are eventually added with other fluxes in
`SumFlux_FluxesOntoGround`. However, I think I'll move them to
waterflux_type because (1) that's needed (or at least consistent with
the above paragraph) if we ever change the state updates to use
individual fluxes rather than summed fluxes, (2) that keeps these
consistent with the mirror-image `qflx_intercepted_rain` and
`qflx_intercepted_snow`, and (3) I can't see a good way to keep these
local without requiring them and/or `qflx_liq_above_canopy` and
`forc_snow_patch` to be stored with a separate tracer dimension. (So I
guess a tentative rule would be: store the more fundamental fluxes in
the types, but can keep various auxiliary / summary things local to the
subroutine that needs them. Here I'm thinking of the `qflx_through_*`
fluxes as being somewhat more "fundamental" than the summed
top-of-canopy inputs.)

One option I could consider in the future is allowing local variables to
have a separate tracer dimension as their last dimension (going from
`bulk_and_tracer_beg` to `bulk_and_tracer_end`). (I may want to check in
with Mariana, and maybe others, about this at some point: okay to do
this, even though it will mean an asymmetry between local variables and
the way we handle variables in water types?) (I thought about introducing
some type to be able to hold a locally-declared water tracer so that the
handling of that would look more like the handling of variables in the
main water types. But I think that may be more complexity than it's
worth.)
