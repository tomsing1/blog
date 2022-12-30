library(paws)
library(purrr)
svc <- dynamodb()

# create the table
svc$create_table(
  AttributeDefinitions = list(
    list(
      AttributeName = "Artist",
      AttributeType = "S"
    ),
    list(
      AttributeName = "SongTitle",
      AttributeType = "S"
    )
  ),
  KeySchema = list(
    list(
      AttributeName = "Artist",
      KeyType = "HASH"
    ),
    list(
      AttributeName = "SongTitle",
      KeyType = "RANGE"
    )
  ),
  ProvisionedThroughput = list(
    ReadCapacityUnits = 5L,
    WriteCapacityUnits = 5L
  ),
  TableName = "Music"
)

# add an item
svc$put_item(
  Item = list(
    AlbumTitle = list(
      S = "Somewhat Famous"
    ),
    Artist = list(
      S = "Acme Band"
    ),
    SongTitle = list(
      S = "Happy Day"
    )
  ),
  ReturnConsumedCapacity = "TOTAL",
  TableName = "Music"
)

# query for the item
response <- svc$query(
  ExpressionAttributeValues = list(
    `:v1` = list(
      S = "Acme Band"
    )
  ),
  KeyConditionExpression = "Artist = :v1",
  ProjectionExpression = "SongTitle",
  TableName = "Music"
)
class(response$Items[[1]]$SongTitle)  # a list
length(response$Items[[1]]$SongTitle)  # one element for each possible data type
purrr::compact(response$Items[[1]]$SongTitle)  # only the `S` field has content

# cleanup
svc$delete_table("Music")