import polars as pl
from google.cloud import bigquery

bq_client = bigquery.Client()


def get_image_logs() -> pl.DataFrame:
    """
    Query bigquery table for image logs.
    Logs are spread across two tables because we had to manually ingest the historical
    logs into a separate table to the one that is written to by the gcp log sink.
    """
    logs_query = """
        select
            timestamp,
            protoPayload.resourceName as full_path,
            protoPayload.authenticationInfo.principalEmail as principal_email,
            protoPayload.authenticationInfo.principalEmail as principal_subject,
            protoPayload.requestMetadata.callerIp as request_ip,
            protoPayload.requestMetadata.callerSuppliedUserAgent as request_user_agent
        from cpg-common.image_logs.historical_logs
        where timestamp < (select min(timestamp) from cpg-common.image_logs.cloudaudit_googleapis_com_data_access)

        UNION ALL

        select
            timestamp,
            protopayload_auditlog.resourceName as full_path,
            protopayload_auditlog.authenticationInfo.principalEmail as principal_email,
            protopayload_auditlog.authenticationInfo.principalSubject as principal_subject,
            protopayload_auditlog.requestMetadata.callerIp as request_ip,
            protopayload_auditlog.requestMetadata.callerSuppliedUserAgent as request_user_agent
        from cpg-common.image_logs.cloudaudit_googleapis_com_data_access
    """
    logs_query_job = bq_client.query(logs_query)
    logs_rows = logs_query_job.result()

    logs_df: pl.DataFrame = pl.from_arrow(logs_rows.to_arrow())
    return (
        # Replace url encoded slashes with real slashes
        logs_df.with_columns(
            pl.col('full_path').str.replace(r'%2F', '/', literal=True),
        )
        # Extract the project, location, repository, image and digest from the image path
        .with_columns(
            pl.col('full_path')
            .str.extract_groups(
                r'projects/(?<image_project>[^/]+)/locations/(?<image_location>[^/]+)/repositories/(?<image_repository>[^/]+)/dockerImages/(?<image_name>[^@]+)@sha256:(?<image_digest>.+)$'
            )
            .alias('path_parts'),
        )
        .unnest('path_parts')
        # remove -archived suffix from repository name so that archived and active
        # images can be compared in logs
        .with_columns(
            # Get the status from the repo name
            pl.when(pl.col('image_repository').str.ends_with('-archive'))
            .then(pl.lit('archived'))
            .otherwise(pl.lit('active'))
            .alias('status'),
            # comptue pre-archive repository
            pl.col('image_repository').str.strip_suffix('-archive'),
            # compute pre-archive short path
            pl.concat_str(
                [
                    pl.col('image_repository').str.strip_suffix('-archive'),
                    pl.col('image_name'),
                ],
                separator='/',
            ).alias('short_path'),
            # compute pre-archive full path
            pl.col('full_path')
            .str.replace('-archive/dockerImages', '/dockerImages', literal=True)
            .alias('full_path'),
        )
    )
