#' Permute membership labels in R (size-preserving).
#' @keywords internal
.mq_perm_labels_R <- function(membership_id) {
    perm <- sample(membership_id, length(membership_id), replace = FALSE)
    names(perm) <- names(membership_id)
    perm
}

#' Update Welford running stats for a vector.
#' @keywords internal
.mq_welford_update <- function(state, x) {
    ok <- is.finite(x) & !is.na(x)
    if (!any(ok)) return(state)
    idx <- which(ok)
    n <- state$n[idx] + 1L
    delta <- x[idx] - state$mean[idx]
    mean <- state$mean[idx] + delta / n
    m2 <- state$m2[idx] + delta * (x[idx] - mean)
    state$n[idx] <- n
    state$mean[idx] <- mean
    state$m2[idx] <- m2
    state
}

#' Finalize Welford stats.
#' @keywords internal
.mq_welford_finalize <- function(state) {
    mean <- state$mean
    sd <- rep(NA_real_, length(mean))
    ok <- state$n > 1L
    sd[ok] <- sqrt(state$m2[ok] / (state$n[ok] - 1L))
    mean[state$n == 0L] <- NA_real_
    list(mean = mean, sd = sd, n = state$n)
}

#' Compute null stats using R label permutations.
#' @keywords internal
.mq_score_null_R <- function(edge_df,
                             membership_id,
                             module_levels,
                             n_perm,
                             seed = 1,
                             compute_conductance = TRUE) {
    if (n_perm <= 0L) {
        return(list(
            member_margin = NULL,
            module_separation = NULL,
            module_conductance = NULL
        ))
    }

    set.seed(seed)
    nodes <- names(membership_id)
    n_nodes <- length(nodes)
    n_modules <- length(module_levels)

    member_state <- list(
        mean = rep(0, n_nodes),
        m2 = rep(0, n_nodes),
        n = rep(0L, n_nodes)
    )
    module_sep_state <- list(
        mean = rep(0, n_modules),
        m2 = rep(0, n_modules),
        n = rep(0L, n_modules)
    )
    module_cond_state <- list(
        mean = rep(0, n_modules),
        m2 = rep(0, n_modules),
        n = rep(0L, n_modules)
    )

    for (i in seq_len(n_perm)) {
        perm <- .mq_perm_labels_R(membership_id)
        member_df <- .mq_member_scores(edge_df, perm, module_levels)
        module_df <- .mq_module_scores(edge_df, perm, module_levels, compute_conductance)

        member_state <- .mq_welford_update(member_state, member_df$member_margin)
        module_sep_state <- .mq_welford_update(module_sep_state, module_df$separation)
        if (isTRUE(compute_conductance)) {
            module_cond_state <- .mq_welford_update(module_cond_state, module_df$conductance)
        }
    }

    list(
        member_margin = .mq_welford_finalize(member_state),
        module_separation = .mq_welford_finalize(module_sep_state),
        module_conductance = if (isTRUE(compute_conductance)) .mq_welford_finalize(module_cond_state) else NULL
    )
}

#' Compute null stats using C++ label permutations if available.
#' @keywords internal
.mq_score_null_cpp <- function(edge_df,
                               membership_id,
                               n_perm,
                               seed = 1,
                               compute_conductance = TRUE) {
    if (!exists("module_quality_perm_null", mode = "function")) {
        return(NULL)
    }
    nodes <- names(membership_id)
    node_index <- seq_along(nodes) - 1L
    names(node_index) <- nodes

    u <- node_index[edge_df$from]
    v <- node_index[edge_df$to]
    w <- edge_df$weight
    keep <- !is.na(u) & !is.na(v) & is.finite(w)
    if (!any(keep)) return(NULL)

    res <- module_quality_perm_null(
        membership = as.integer(membership_id[nodes]),
        u = as.integer(u[keep]),
        v = as.integer(v[keep]),
        w = as.numeric(w[keep]),
        n_perm = as.integer(n_perm),
        seed = as.integer(seed),
        do_member = TRUE,
        do_module = TRUE,
        do_conductance = isTRUE(compute_conductance)
    )

    list(
        member_margin = list(mean = res$member_margin_mean, sd = res$member_margin_sd),
        module_separation = list(mean = res$module_sep_mean, sd = res$module_sep_sd),
        module_conductance = if (isTRUE(compute_conductance)) {
            list(mean = res$module_cond_mean, sd = res$module_cond_sd)
        } else {
            NULL
        }
    )
}

#' Compute z/p/fdr from null mean/sd.
#' @keywords internal
.mq_z_p_from_null <- function(observed, null_mean, null_sd) {
    z <- rep(NA_real_, length(observed))
    p <- rep(NA_real_, length(observed))
    ok <- is.finite(observed) & is.finite(null_mean) & is.finite(null_sd) & null_sd > 0
    if (any(ok)) {
        z[ok] <- (observed[ok] - null_mean[ok]) / null_sd[ok]
        p[ok] <- 2 * stats::pnorm(abs(z[ok]), lower.tail = FALSE)
    }
    fdr <- stats::p.adjust(p, method = "fdr")
    list(z = z, p = p, fdr = fdr)
}
