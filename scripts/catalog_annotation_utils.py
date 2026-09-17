"""Helpers shared by the scripts that add per-locus annotations to a TR catalog JSON file.

Each of those scripts owns a fixed set of annotation fields and has to replace them on every run:
a locus that drops out of the script's input table must lose the values an earlier run gave it,
rather than keeping them forever. These helpers keep that behavior the same across scripts.
"""


def clear_previous_annotations(record, fields):
    """Remove the calling script's own annotation fields from a catalog record.

    Args:
        record: catalog record dictionary. Modified in place.
        fields: the annotation field names the calling script owns.

    Returns:
        True if the record carried any of those fields from an earlier run.
    """
    # Pop every field before testing any of them. A generator inside any() would short-circuit on
    # the first field that was present and leave the rest of a stale annotation behind.
    previous_values = [record.pop(field, None) for field in fields]
    return any(value is not None for value in previous_values)


def find_unlisted_annotation_fields(annotation_lookup, fields):
    """Return the annotation fields that are about to be written but are missing from `fields`.

    A field the calling script doesn't list is never cleared, so it would go stale on a re-run.

    Args:
        annotation_lookup: dictionary of locus id -> the annotations to write for that locus.
        fields: the annotation field names the calling script owns.

    Returns:
        The set of field names that are written but unlisted. Empty when the list is complete.
    """
    return {field for annotations in annotation_lookup.values() for field in annotations} - set(fields)


def print_annotation_replacement_summary(label, loss_reason, cleared_locus_count, lost_locus_count):
    """Report how many loci had earlier annotations replaced, and how many lost them entirely.

    Both lines stay silent on a normal pipeline run, where the catalog arrives unannotated.

    Args:
        label: name of the annotation family, e.g. "AoU1027" or "variation cluster".
        loss_reason: why the lost loci got nothing back, phrased to follow "because they",
            e.g. "are no longer in /path/to/table.tsv.gz". Leaving a locus unannotated isn't
            always the same as dropping it from the input table, so the caller words this.
        cleared_locus_count: loci that carried these annotations from an earlier run.
        lost_locus_count: how many of those got nothing back.
    """
    if cleared_locus_count > 0:
        print(f"  - {cleared_locus_count:,d} already carried {label} annotations from an earlier run, "
              f"which were discarded before re-annotating")
    if lost_locus_count > 0:
        print(f"  - {lost_locus_count:,d} of those lost their {label} annotations entirely because "
              f"they {loss_reason}")
