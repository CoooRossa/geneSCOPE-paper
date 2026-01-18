#' Finalize member-level output table.
#' @keywords internal
.mq_make_run_id <- function() {
    paste0(format(Sys.time(), "%Y%m%d%H%M%S"), "_", sprintf("%05d", Sys.getpid()))
}

#' Finalize member-level output table.
#' @keywords internal
.mq_finalize_member_table <- function(member_scores) {
    if (is.null(member_scores) || !nrow(member_scores)) {
        return(data.frame())
    }
    member_scores$score_member <- member_scores$member_margin
    member_scores$z_member <- member_scores$member_z
    member_scores$p_member <- member_scores$member_p
    member_scores$padj_member <- member_scores$member_fdr
    member_scores$n_intra_edges <- member_scores$within_n
    member_scores$n_inter_edges <- member_scores$between_n
    member_scores$sum_intra_w <- member_scores$within_sum
    member_scores$sum_inter_w <- member_scores$between_sum

    rank_fun <- function(x) {
        rank(-x, ties.method = "min", na.last = "keep")
    }
    member_scores$rank_in_module <- ave(member_scores$score_member, member_scores$module, FUN = rank_fun)
    member_scores
}

#' Build top-members summary per module.
#' @keywords internal
.mq_top_members <- function(member_scores, top_n = 5L) {
    if (is.null(member_scores) || !nrow(member_scores)) {
        return(setNames(character(0), character(0)))
    }
    top_n <- max(1L, as.integer(top_n))
    by_module <- split(member_scores, member_scores$module)
    top_members <- vapply(by_module, function(df) {
        ord <- order(df$score_member, decreasing = TRUE, na.last = TRUE)
        genes <- df$gene[ord]
        genes <- genes[seq_len(min(top_n, length(genes)))]
        paste(genes, collapse = ",")
    }, character(1))
    top_members
}

#' Finalize module-level output table.
#' @keywords internal
.mq_finalize_module_table <- function(module_scores,
                                      member_scores,
                                      null_summary = NULL,
                                      top_n = 5L) {
    if (is.null(module_scores) || !nrow(module_scores)) {
        return(data.frame())
    }
    module_scores$n_members <- module_scores$module_size
    module_scores$score_module <- module_scores$separation
    module_scores$effect_size <- module_scores$separation
    module_scores$p_module <- module_scores$separation_p
    module_scores$padj_module <- module_scores$separation_fdr

    if (!("null_mean" %in% names(module_scores)) || !("null_sd" %in% names(module_scores))) {
        if (!is.null(null_summary) && !is.null(null_summary$module_separation)) {
            module_scores$null_mean <- null_summary$module_separation$mean
            module_scores$null_sd <- null_summary$module_separation$sd
        } else {
            module_scores$null_mean <- NA_real_
            module_scores$null_sd <- NA_real_
        }
    }

    top_members <- .mq_top_members(member_scores, top_n = top_n)
    module_scores$top_members <- top_members[as.character(module_scores$module)]
    module_scores
}

#' Create summary tables (module + member top-N).
#' @keywords internal
.mq_summary_tables <- function(member_scores, module_scores, top_n = 5L) {
    member_scores <- .mq_finalize_member_table(member_scores)
    module_scores <- .mq_finalize_module_table(module_scores, member_scores, null_summary = NULL, top_n = top_n)

    member_table <- data.frame()
    if (nrow(member_scores)) {
        member_table <- do.call(rbind, lapply(split(member_scores, member_scores$module), function(df) {
            ord <- order(df$score_member, decreasing = TRUE, na.last = TRUE)
            df <- df[ord, , drop = FALSE]
            df[seq_len(min(top_n, nrow(df))), , drop = FALSE]
        }))
    }

    module_table <- module_scores
    if (nrow(module_table)) {
        ord <- order(module_table$padj_module, -module_table$score_module, na.last = TRUE)
        module_table <- module_table[ord, , drop = FALSE]
    }

    list(module_table = module_table, member_table = member_table)
}
